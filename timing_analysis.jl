# Focused timing analysis to identify CLBM bottlenecks

include("src/CLBM/clbm_config.jl")

using SymPy
using LinearAlgebra
using SparseArrays

# Include files one by one to time each
println("=== CLBM TIMING ANALYSIS ===")
println("Configuration: ngrid=$ngrid, use_sparse=$use_sparse")
println()

# Time individual file includes to see setup overhead
print("Loading coeffs_poly... ")
@time include("src/CLBM/coeffs_poly.jl")

print("Loading collision_sym... ")
@time include("src/CLBM/collision_sym.jl")

print("Loading carleman_transferA... ")
@time include("src/CLBM/carleman_transferA.jl")

print("Loading carleman_transferA_ngrid... ")
@time include("src/CLBM/carleman_transferA_ngrid.jl")

print("Loading LBM_const_subs... ")
@time include("src/CLBM/LBM_const_subs.jl")

print("Loading lbm_cons... ")
@time include("src/LBM/lbm_cons.jl")

print("Loading lbm_const_sym... ")
@time include("src/LBM/lbm_const_sym.jl")

print("Loading forcing... ")
@time include("src/LBM/forcing.jl")

print("Loading f_initial... ")
@time include("src/LBM/f_initial.jl")

print("Loading timeMarching... ")
@time include("src/CLBM/timeMarching.jl")

print("Loading CLBM_collision_test... ")
@time include("src/CLBM/CLBM_collision_test.jl")

println("\n=== FUNCTION TIMING ===")

# Time individual function calls
print("lbm_const_sym()... ")
const_time = @elapsed begin
    w, e, w_val, e_val = lbm_const_sym()
    global w_value = w_val
    global e_value = e_val
end
println("$(round(const_time, digits=3))s")

print("collision()... ")
collision_time = @elapsed begin
    f, omega, u, rho = collision(Q, D, w, e, rho0, lTaylor, lorder2)
end
println("$(round(collision_time, digits=3))s")

print("get_coeff_LBM_Fi_ngrid()... ")
coeff_time = @elapsed begin
    global F1_ngrid, F2_ngrid, F3_ngrid = get_coeff_LBM_Fi_ngrid(poly_order, Q, f, omega, tau_value, ngrid)
end
println("$(round(coeff_time, digits=3))s")

# Matrix construction timing
if use_sparse
    print("carleman_C_sparse()... ")
    matrix_time = @elapsed begin
        C_sparse, bt_sparse, F0_sparse = carleman_C_sparse(Q, truncation_order, poly_order, f, omega, tau_value, force_factor, w_value, e_value)
    end
    println("$(round(matrix_time, digits=3))s")
    
    print("Matrix properties: ")
    println("$(size(C_sparse, 1))×$(size(C_sparse, 2)), nnz=$(nnz(C_sparse)), sparsity=$(round((1-nnz(C_sparse)/prod(size(C_sparse)))*100, digits=1))%")
else
    print("carleman_C()... ")
    matrix_time = @elapsed begin
        C_dense, bt_dense, F0_dense = carleman_C(Q, truncation_order, poly_order, f, omega, tau_value, force_factor, w_value, e_value)
    end
    println("$(round(matrix_time, digits=3))s")
    
    print("Matrix properties: ")
    println("$(size(C_dense, 1))×$(size(C_dense, 2))")
end

# Time stepping analysis
println("\n=== TIME STEPPING ANALYSIS ===")

# Test with different time step counts
test_steps = [5, 10, 25]

for steps in test_steps
    print("$steps time steps... ")
    
    if use_sparse
        step_time = @elapsed begin
            fT, VT_f, VT = CLBM_collision_test_sparse(Q, omega, f, truncation_order, dt, tau_value, e_value, steps, false)
        end
    else
        step_time = @elapsed begin
            fT, VT_f, VT = CLBM_collision_test(Q, omega, f, C_dense, truncation_order, dt, tau_value, e_value, steps, false)
        end
    end
    
    ms_per_step = step_time * 1000 / steps
    println("$(round(step_time, digits=3))s ($(round(ms_per_step, digits=2)) ms/step)")
end

# Detailed breakdown of major components
println("\n=== COMPONENT BREAKDOWN ===")
total_setup = const_time + collision_time + coeff_time + matrix_time
println("Setup phase breakdown:")
println("  LBM constants: $(round(const_time/total_setup*100, digits=1))%")
println("  Collision operators: $(round(collision_time/total_setup*100, digits=1))%")
println("  Coefficient computation: $(round(coeff_time/total_setup*100, digits=1))%")
println("  Matrix construction: $(round(matrix_time/total_setup*100, digits=1))%")

println("\nTotal setup time: $(round(total_setup, digits=3))s")

# Memory usage estimate
if use_sparse
    matrix_memory = (length(C_sparse.nzval) * 8 + length(C_sparse.rowval) * 8 + length(C_sparse.colptr) * 8) / 1024^2
    println("Sparse matrix memory: $(round(matrix_memory, digits=2)) MB")
else
    matrix_memory = sizeof(C_dense) / 1024^2
    println("Dense matrix memory: $(round(matrix_memory, digits=2)) MB")
end

println("\n=== PERFORMANCE RECOMMENDATIONS ===")
if matrix_time > const_time + collision_time + coeff_time
    println("🔍 BOTTLENECK: Matrix construction ($(round(matrix_time, digits=3))s)")
    if use_sparse
        println("   - Consider optimizing sparse assembly algorithms")
        println("   - Check if block overlap handling can be optimized")
    else
        println("   - Consider switching to sparse matrices for ngrid > 1")
        println("   - Dense construction scales poorly with problem size")
    end
elseif collision_time > matrix_time
    println("🔍 BOTTLENECK: Collision operator computation ($(round(collision_time, digits=3))s)")
    println("   - Consider caching symbolic computations")
    println("   - Check if SymPy operations can be optimized")
elseif coeff_time > matrix_time
    println("🔍 BOTTLENECK: Coefficient computation ($(round(coeff_time, digits=3))s)")
    println("   - Consider pre-computing F matrices")
    println("   - Check get_coeff_LBM_Fi_ngrid efficiency")
else
    println("✅ No major bottlenecks in setup phase")
    println("   Time stepping performance: $(round(ms_per_step, digits=2)) ms/step")
end

if use_sparse && matrix_memory > 100
    println("⚠️  Large sparse matrix memory usage: $(round(matrix_memory, digits=1)) MB")
    println("   Consider further sparsity optimizations")
end

println("\n=== TIMING COMPLETE ===")
println("Run this script to identify where time is spent in your CLBM simulation.")
