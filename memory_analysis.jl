# Memory analysis for different ngrid values

include("src/CLBM/clbm_config.jl")

using LinearAlgebra

# Load just the dimension calculation function
include("src/CLBM/carleman_transferA.jl")

println("=== MEMORY ANALYSIS FOR DIFFERENT NGRID VALUES ===")

function analyze_memory_requirements(ngrid_val, Q=3, truncation_order=3)
    # Calculate Carleman matrix dimension
    C_dim = carleman_C_dim(Q, truncation_order, ngrid_val)
    
    # Calculate memory requirements
    total_elements = C_dim^2
    dense_memory_bytes = total_elements * 8  # 8 bytes per Float64
    dense_memory_mb = dense_memory_bytes / 1024^2
    dense_memory_gb = dense_memory_mb / 1024
    
    # Estimate sparse matrix memory (very rough estimate)
    # For Carleman matrices, sparsity typically increases with problem size
    if ngrid_val == 1
        sparsity_estimate = 0.8  # 80% sparse
    elseif ngrid_val == 2
        sparsity_estimate = 0.999  # 99.9% sparse (observed)
    else
        # Extrapolate - larger problems tend to be sparser
        sparsity_estimate = min(0.9999, 0.999 + (ngrid_val - 2) * 0.0001)
    end
    
    nnz_estimate = Int(ceil(total_elements * (1 - sparsity_estimate)))
    sparse_memory_bytes = nnz_estimate * 8 + nnz_estimate * 8 + C_dim * 8  # values + indices + colptr
    sparse_memory_mb = sparse_memory_bytes / 1024^2
    
    return (
        ngrid = ngrid_val,
        matrix_size = C_dim,
        total_elements = total_elements,
        dense_memory_mb = dense_memory_mb,
        dense_memory_gb = dense_memory_gb,
        estimated_sparsity = sparsity_estimate * 100,
        estimated_nnz = nnz_estimate,
        sparse_memory_mb = sparse_memory_mb
    )
end

# Test different ngrid values
test_values = [1, 2, 3, 4, 5, 6, 7, 8, 10]

println("ngrid | Matrix Size | Dense Memory | Sparse Memory | Sparsity")
println("------|-------------|--------------|---------------|----------")

memory_critical = []
memory_impossible = []

for ngrid_val in test_values
    try
        result = analyze_memory_requirements(ngrid_val)
        
        # Format output
        if result.dense_memory_gb < 1
            dense_str = "$(round(result.dense_memory_mb, digits=1)) MB"
        else
            dense_str = "$(round(result.dense_memory_gb, digits=1)) GB"
        end
        
        sparse_str = "$(round(result.sparse_memory_mb, digits=1)) MB"
        sparsity_str = "$(round(result.estimated_sparsity, digits=1))%"
        
        println("$(lpad(ngrid_val, 5)) | $(lpad(result.matrix_size, 11)) | $(lpad(dense_str, 12)) | $(lpad(sparse_str, 13)) | $(lpad(sparsity_str, 8))")
        
        # Flag problematic cases
        if result.dense_memory_gb > 8  # More than 8GB
            push!(memory_impossible, ngrid_val)
        elseif result.dense_memory_gb > 1  # More than 1GB
            push!(memory_critical, ngrid_val)
        end
        
    catch e
        println("$(lpad(ngrid_val, 5)) | ERROR: $e")
        push!(memory_impossible, ngrid_val)
    end
end

println("\n=== MEMORY ANALYSIS SUMMARY ===")

if !isempty(memory_critical)
    println("⚠️  MEMORY CRITICAL (>1GB dense): ngrid = $(join(memory_critical, ", "))")
    println("   These values require sparse matrices to be feasible")
end

if !isempty(memory_impossible)
    println("❌ MEMORY IMPOSSIBLE (>8GB dense): ngrid = $(join(memory_impossible, ", "))")
    println("   These values may cause out-of-memory even with sparse matrices")
end

# Specific analysis for ngrid=7
println("\n=== DETAILED ANALYSIS FOR NGRID=7 ===")
try
    result = analyze_memory_requirements(7)
    println("Matrix dimension: $(result.matrix_size) × $(result.matrix_size)")
    println("Total matrix elements: $(result.total_elements)")
    println("Dense memory required: $(round(result.dense_memory_gb, digits=1)) GB")
    println("Estimated sparse memory: $(round(result.sparse_memory_mb, digits=1)) MB")
    println("Estimated sparsity: $(round(result.estimated_sparsity, digits=2))%")
    
    if result.dense_memory_gb > 8
        println("\n🚨 DIAGNOSIS: ngrid=7 requires $(round(result.dense_memory_gb, digits=1)) GB of RAM")
        println("   This exceeds typical system memory, causing the OS to kill the process")
        println("   Even sparse matrices may be challenging at this size")
    end
    
    # Check current configuration
    println("\nCurrent configuration check:")
    if use_sparse
        println("✅ use_sparse = true (good for ngrid=7)")
        println("   But the matrix may still be too large even in sparse format")
    else
        println("❌ use_sparse = false (will definitely fail for ngrid=7)")
        println("   Dense matrices are impossible at this size")
    end
    
catch e
    println("ERROR analyzing ngrid=7: $e")
end

println("\n=== RECOMMENDATIONS ===")
println("1. For ngrid ≥ 3: ALWAYS use sparse matrices (use_sparse = true)")
println("2. For ngrid ≥ 5: Consider reducing truncation_order if possible")
println("3. For ngrid ≥ 7: May need algorithmic changes or more RAM")
println("4. Monitor system memory usage with 'htop' or Activity Monitor")

println("\n=== SAFE TESTING PROCEDURE ===")
println("To test large ngrid values safely:")
println("1. Ensure use_sparse = true")
println("2. Start with smaller truncation_order (e.g., 2 instead of 3)")
println("3. Monitor memory usage during matrix construction")
println("4. Test matrix construction before time stepping:")
println("   julia -e 'include(\"timing_analysis.jl\")' | head -20")

# Check system memory if possible
try
    # Try to get system memory info (works on some systems)
    if Sys.islinux()
        memory_info = read("/proc/meminfo", String)
        if occursin("MemTotal", memory_info)
            total_kb = parse(Int, match(r"MemTotal:\s*(\d+)", memory_info).captures[1])
            total_gb = total_kb / 1024^2
            println("\nSystem memory: ~$(round(total_gb, digits=1)) GB")
        end
    elseif Sys.isapple()
        println("\nTo check system memory on macOS: 'sysctl hw.memsize'")
    end
catch
    println("\nTo check available system memory:")
    println("  macOS: sysctl hw.memsize")
    println("  Linux: cat /proc/meminfo | grep MemTotal")
end
