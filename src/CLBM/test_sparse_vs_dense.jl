# Test to verify sparse and dense Carleman matrix implementations give identical results

l_sympy = true
QCFD_SRC = ENV["QCFD_SRC"]  
QCFD_HOME = ENV["QCFD_HOME"]  

using Test
using LinearAlgebra

# Load centralized configuration
include("clbm_config.jl")

include(QCFD_HOME * "/visualization/plot_kit.jl")

if l_sympy
    using SymPy
    using LinearAlgebra
    include(QCFD_SRC * "CLBM/coeffs_poly.jl")
else
    using Symbolics
end

# Include necessary files
include(QCFD_SRC * "CLBM/collision_sym.jl")
include(QCFD_SRC * "CLBM/carleman_transferA.jl")
include(QCFD_SRC * "CLBM/carleman_transferA_ngrid.jl")
include(QCFD_SRC * "CLBM/LBM_const_subs.jl")
include(QCFD_SRC * "LBM/lbm_cons.jl")
include(QCFD_SRC * "LBM/lbm_const_sym.jl")
include(QCFD_SRC * "LBM/forcing.jl")
include(QCFD_SRC * "LBM/f_initial.jl")
include(QCFD_SRC * "CLBM/timeMarching.jl")

function test_sparse_vs_dense_carleman()
    println("Testing sparse vs dense Carleman matrix implementations...")
    
    # Override n_time for quick testing
    local_n_time = 10  # Small number for quick test
    
    # Set up LBM constants (updates global w_value, e_value)
    w, e, w_val, e_val = lbm_const_sym()
    global w_value = w_val
    global e_value = e_val
    
    # Generate collision operators
    f, omega, u, rho = collision(Q, D, w, e, rho0, lTaylor, lorder2)
    
    # Set up initial conditions
    f_ini = f_ini_test(u0)
    
    # Initialize global variables needed for transferA_ngrid
    global F1_ngrid, F2_ngrid, F3_ngrid = get_coeff_LBM_Fi_ngrid(poly_order, Q, f, omega, tau_value, ngrid)
    
    # Get the dense Carleman matrix
    C_dense, bt_dense, F0_dense = carleman_C(Q, truncation_order, poly_order, f, omega, tau_value, force_factor, w_val, e_val)
    
    # Get the sparse Carleman matrix
    C_sparse, bt_sparse, F0_sparse = carleman_C_sparse(Q, truncation_order, poly_order, f, omega, tau_value, force_factor, w_val, e_val)
    
    # Test 1: Compare matrices
    println("Test 1: Comparing Carleman matrices...")
    C_sparse_full = Array(C_sparse)  # Convert sparse to dense for comparison
    @test isapprox(C_dense, C_sparse_full, rtol=1e-12)
    println("✓ Carleman matrices match within tolerance")
    
    # Test 2: Compare forcing vectors
    println("Test 2: Comparing forcing vectors...")
    bt_sparse_full = Array(bt_sparse)  # Convert sparse to dense for comparison
    @test isapprox(bt_dense, bt_sparse_full, rtol=1e-12)
    @test isapprox(F0_dense, F0_sparse, rtol=1e-12)
    println("✓ Forcing vectors match within tolerance")
    
    # Test 3: Compare time marching results
    println("Test 3: Comparing time marching results...")
    
    # Run dense version
    VT_f_dense, VT_dense, uT_dense, fT_dense = timeMarching_collision_CLBM(
        omega, f, tau_value, Q, C_dense, truncation_order, e_val, dt, f_ini, local_n_time, false
    )
    
    # Run sparse version
    VT_f_sparse, VT_sparse, uT_sparse, fT_sparse = timeMarching_collision_CLBM_sparse(
        omega, f, tau_value, Q, truncation_order, e_val, dt, f_ini, local_n_time, false
    )
    
    # Compare results
    @test isapprox(VT_f_dense, VT_f_sparse, rtol=1e-10)
    @test isapprox(VT_dense, VT_sparse, rtol=1e-10)
    @test isapprox(uT_dense, uT_sparse, rtol=1e-10)
    @test isapprox(fT_dense, fT_sparse, rtol=1e-10)
    
    println("✓ Time marching results match within tolerance")
    
    # Test 4: Performance comparison
    println("Test 4: Performance comparison...")
    
    # Time dense version
    time_dense = @elapsed begin
        for i = 1:5
            timeMarching_collision_CLBM(
                omega, f, tau_value, Q, C_dense, truncation_order, e_val, dt, f_ini, local_n_time, false
            )
        end
    end
    
    # Time sparse version
    time_sparse = @elapsed begin
        for i = 1:5
            timeMarching_collision_CLBM_sparse(
                omega, f, tau_value, Q, truncation_order, e_val, dt, f_ini, local_n_time, false
            )
        end
    end
    
    println("Dense version average time: $(round(time_dense/5, digits=4)) seconds")
    println("Sparse version average time: $(round(time_sparse/5, digits=4)) seconds")
    
    if time_sparse < time_dense
        println("✓ Sparse version is faster!")
        speedup = time_dense / time_sparse
        println("Speedup: $(round(speedup, digits=2))x")
    else
        println("Note: For small problems, sparse overhead may make it slower")
    end
    
    # Test 5: Memory usage comparison
    println("Test 5: Memory usage comparison...")
    
    C_dense_memory = sizeof(C_dense) / 1024^2  # MB
    C_sparse_memory = (length(C_sparse.nzval) * sizeof(Float64) + 
                       length(C_sparse.rowval) * sizeof(Int) + 
                       length(C_sparse.colptr) * sizeof(Int)) / 1024^2  # MB
    
    println("Dense matrix memory: $(round(C_dense_memory, digits=2)) MB")
    println("Sparse matrix memory: $(round(C_sparse_memory, digits=2)) MB")
    
    memory_savings = (C_dense_memory - C_sparse_memory) / C_dense_memory * 100
    println("Memory savings: $(round(memory_savings, digits=1))%")
    
    # Test 6: Sparsity analysis
    println("Test 6: Sparsity analysis...")
    
    total_elements = size(C_dense, 1) * size(C_dense, 2)
    nonzero_elements = nnz(C_sparse)
    sparsity = (total_elements - nonzero_elements) / total_elements * 100
    
    println("Matrix size: $(size(C_dense, 1)) × $(size(C_dense, 2))")
    println("Total elements: $total_elements")
    println("Non-zero elements: $nonzero_elements")
    println("Sparsity: $(round(sparsity, digits=1))%")
    
    println("\n🎉 All tests passed! Sparse and dense implementations are equivalent.")
    
    return true
end

function analyze_sparsity_metrics(C_dense, C_sparse, ngrid_val)
    """Extract sparsity and memory metrics from matrices"""
    total_elements = size(C_dense, 1) * size(C_dense, 2)
    nonzero_elements = nnz(C_sparse)
    sparsity = (total_elements - nonzero_elements) / total_elements * 100
    
    dense_memory = sizeof(C_dense) / 1024^2  # MB
    sparse_memory = (length(C_sparse.nzval) * sizeof(Float64) + 
                     length(C_sparse.rowval) * sizeof(Int) + 
                     length(C_sparse.colptr) * sizeof(Int)) / 1024^2  # MB
    
    memory_savings = (dense_memory - sparse_memory) / dense_memory * 100
    
    return Dict(
        "ngrid" => ngrid_val,
        "matrix_size" => size(C_dense),
        "total_elements" => total_elements,
        "nonzero_elements" => nonzero_elements,
        "sparsity" => sparsity,
        "dense_memory" => dense_memory,
        "sparse_memory" => sparse_memory,
        "memory_savings" => memory_savings
    )
end

function test_ngrid_sparsity_comparison()
    println("\n" * "="^60)
    println("TESTING SPARSITY IMPROVEMENTS WITH GRID SIZE")
    println("="^60)
    
    # Store original ngrid
    original_ngrid = ngrid
    results = []
    
    for test_ngrid in [1, 2]
        println("\nTesting ngrid = $test_ngrid...")
        
        # Update global ngrid
        global ngrid = test_ngrid
        
        # Set up LBM constants
        w, e, w_val, e_val = lbm_const_sym()
        global w_value = w_val
        global e_value = e_val
        
        # Generate collision operators
        f, omega, u, rho = collision(Q, D, w, e, rho0, lTaylor, lorder2)
        
        # Update ngrid-dependent coefficients
        global F1_ngrid, F2_ngrid, F3_ngrid = get_coeff_LBM_Fi_ngrid(poly_order, Q, f, omega, tau_value, test_ngrid)
        
        # Get matrices
        C_dense, _, _ = carleman_C(Q, truncation_order, poly_order, f, omega, tau_value, force_factor, w_val, e_val)
        C_sparse, _, _ = carleman_C_sparse(Q, truncation_order, poly_order, f, omega, tau_value, force_factor, w_val, e_val)
        
        # Analyze metrics
        metrics = analyze_sparsity_metrics(C_dense, C_sparse, test_ngrid)
        push!(results, metrics)
        
        println("  Matrix size: $(metrics["matrix_size"][1]) × $(metrics["matrix_size"][2])")
        println("  Sparsity: $(round(metrics["sparsity"], digits=1))%")
        println("  Memory savings: $(round(metrics["memory_savings"], digits=1))%")
    end
    
    # Restore original ngrid
    global ngrid = original_ngrid
    
    # Compare results
    println("\n" * "="^40)
    println("COMPARISON RESULTS")
    println("="^40)
    
    ngrid1_result = results[1]
    ngrid2_result = results[2]
    
    size_growth = ngrid2_result["total_elements"] / ngrid1_result["total_elements"]
    sparsity_increase = ngrid2_result["sparsity"] - ngrid1_result["sparsity"]
    savings_improvement = ngrid2_result["memory_savings"] - ngrid1_result["memory_savings"]
    
    println("Matrix size growth (ngrid=2 vs ngrid=1): $(round(size_growth, digits=2))x")
    println("Sparsity increase: +$(round(sparsity_increase, digits=1)) percentage points")
    println("Memory savings improvement: +$(round(savings_improvement, digits=1)) percentage points")
    
    # Validation
    if sparsity_increase > 0
        println("✓ CONFIRMED: Sparse matrices become more beneficial for larger grids")
    else
        println("⚠ WARNING: Expected sparsity increase not observed")
    end
    
    return results
end

# Run the tests
if abspath(PROGRAM_FILE) == @__FILE__
    test_sparse_vs_dense_carleman()
    test_ngrid_sparsity_comparison()
end
