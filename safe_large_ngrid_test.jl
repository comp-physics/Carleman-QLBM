# Safe testing script for large ngrid values

include("src/CLBM/clbm_config.jl")

println("=== SAFE LARGE NGRID TEST ===")
println("Testing ngrid = $ngrid")

# Force sparse mode for safety
if !use_sparse
    println("⚠️  WARNING: Forcing use_sparse = true for ngrid = $ngrid")
    global use_sparse = true
end

using LinearAlgebra

# Load minimal dependencies first
include("src/CLBM/carleman_transferA.jl")

# Calculate expected memory usage
function estimate_memory_gb(ngrid_val, Q=3, truncation_order=3)
    C_dim = carleman_C_dim(Q, truncation_order, ngrid_val)
    total_elements = C_dim^2
    dense_memory_gb = total_elements * 8 / 1024^3
    return C_dim, dense_memory_gb
end

matrix_size, estimated_gb = estimate_memory_gb(ngrid)

println("Matrix size will be: $matrix_size × $matrix_size")
println("Dense memory would require: $(round(estimated_gb, digits=1)) GB")

# Safety checks
if estimated_gb > 16
    println("🚨 ABORT: Matrix too large (>16GB)")
    println("   This will likely cause system to freeze or kill the process")
    println("   Consider reducing ngrid or truncation_order")
    exit(1)
elseif estimated_gb > 4
    println("⚠️  WARNING: Large matrix ($(round(estimated_gb, digits=1)) GB)")
    println("   Proceeding with caution using sparse matrices only...")
    
    # Ask for confirmation (in interactive mode)
    if isinteractive()
        print("Continue? (y/N): ")
        response = readline()
        if lowercase(strip(response)) != "y"
            println("Aborted by user")
            exit(0)
        end
    else
        println("   Non-interactive mode: proceeding automatically")
    end
end

# Load remaining dependencies
println("\nLoading remaining dependencies...")

using SymPy
using SparseArrays

include("src/CLBM/coeffs_poly.jl")
include("src/CLBM/collision_sym.jl")
include("src/CLBM/carleman_transferA_ngrid.jl")
include("src/CLBM/LBM_const_subs.jl")
include("src/LBM/lbm_cons.jl")
include("src/LBM/lbm_const_sym.jl")
include("src/LBM/forcing.jl")
include("src/LBM/f_initial.jl")
include("src/CLBM/timeMarching.jl")
include("src/CLBM/CLBM_collision_test.jl")

println("✅ Dependencies loaded")

# Phase 1: Setup with timing
println("\n=== PHASE 1: SETUP ===")

print("Setting up LBM constants... ")
@time begin
    w, e, w_val, e_val = lbm_const_sym()
    global w_value = w_val
    global e_value = e_val
end

print("Computing collision operators... ")
@time begin
    f, omega, u, rho = collision(Q, D, w, e, rho0, lTaylor, lorder2)
end

print("Computing F coefficients... ")
@time begin
    global F1_ngrid, F2_ngrid, F3_ngrid = get_coeff_LBM_Fi_ngrid(poly_order, Q, f, omega, tau_value, ngrid)
end

println("✅ Setup complete")

# Phase 2: Matrix construction with memory monitoring
println("\n=== PHASE 2: MATRIX CONSTRUCTION ===")
println("Attempting sparse matrix construction...")

# Try to get memory usage before matrix construction
function get_memory_usage_mb()
    try
        if Sys.islinux()
            # Parse /proc/meminfo for available memory
            meminfo = read("/proc/meminfo", String)
            available_match = match(r"MemAvailable:\s*(\d+)", meminfo)
            if available_match !== nothing
                return parse(Int, available_match.captures[1]) / 1024  # Convert KB to MB
            end
        end
        return nothing
    catch
        return nothing
    end
end

initial_memory = get_memory_usage_mb()
if initial_memory !== nothing
    println("Available memory before construction: $(round(initial_memory/1024, digits=1)) GB")
end

# Construct matrix with error handling
matrix_construction_time = 0.0
success = false

try
    print("Building sparse matrix... ")
    matrix_construction_time = @elapsed begin
        C_sparse, bt_sparse, F0_sparse = carleman_C_sparse(Q, truncation_order, poly_order, f, omega, tau_value, force_factor, w_value, e_value)
    end
    
    # Check matrix properties
    actual_size = size(C_sparse)
    actual_nnz = nnz(C_sparse)
    actual_sparsity = (1 - actual_nnz / prod(actual_size)) * 100
    actual_memory_mb = (length(C_sparse.nzval) * 8 + length(C_sparse.rowval) * 8 + length(C_sparse.colptr) * 8) / 1024^2
    
    println("$(round(matrix_construction_time, digits=2))s")
    println("✅ Matrix construction successful!")
    println("   Actual size: $(actual_size[1]) × $(actual_size[2])")
    println("   Non-zeros: $actual_nnz")
    println("   Sparsity: $(round(actual_sparsity, digits=2))%")
    println("   Memory usage: $(round(actual_memory_mb, digits=1)) MB")
    
    success = true
    
catch e
    println("FAILED")
    println("❌ Matrix construction failed: $e")
    
    if isa(e, OutOfMemoryError) || occursin("memory", string(e))
        println("   This appears to be a memory-related failure")
        println("   Try reducing ngrid or truncation_order")
    end
    
    exit(1)
end

final_memory = get_memory_usage_mb()
if final_memory !== nothing && initial_memory !== nothing
    memory_used = initial_memory - final_memory
    println("Estimated memory consumed: $(round(memory_used, digits=1)) MB")
end

# Phase 3: Short time stepping test
if success
    println("\n=== PHASE 3: TIME STEPPING TEST ===")
    
    # Test with very few time steps first
    test_steps = 3
    println("Testing $test_steps time steps...")
    
    try
        step_time = @elapsed begin
            fT, VT_f, VT = CLBM_collision_test_sparse(Q, omega, f, truncation_order, dt, tau_value, e_value, test_steps, false)
        end
        
        println("✅ Time stepping successful!")
        println("   Time for $test_steps steps: $(round(step_time, digits=3))s")
        println("   Time per step: $(round(step_time*1000/test_steps, digits=1)) ms")
        
        # Test scaling with more steps if the first test was fast enough
        if step_time < 1.0  # If 3 steps took less than 1 second
            println("\nTesting longer simulation...")
            longer_steps = 10
            longer_time = @elapsed begin
                fT, VT_f, VT = CLBM_collision_test_sparse(Q, omega, f, truncation_order, dt, tau_value, e_value, longer_steps, false)
            end
            println("✅ $longer_steps steps completed in $(round(longer_time, digits=3))s")
            println("   Projected time for n_time=$n_time: $(round(longer_time * n_time / longer_steps, digits=1))s")
        else
            println("⚠️  Time stepping is slow for this matrix size")
        end
        
    catch e
        println("❌ Time stepping failed: $e")
    end
end

println("\n=== TEST COMPLETE ===")
println("Matrix construction time: $(round(matrix_construction_time, digits=2))s")

if success
    println("✅ ngrid=$ngrid is feasible with sparse matrices")
    println("   You can now run the full simulation with:")
    println("   julia src/CLBM/clbm_run.jl")
else
    println("❌ ngrid=$ngrid is not feasible with current resources")
    println("   Consider reducing the problem size")
end
