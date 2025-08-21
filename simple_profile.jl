# Simple profiling to identify CLBM bottlenecks

include("src/CLBM/clbm_config.jl")

using SymPy
using LinearAlgebra
using SparseArrays

include("src/CLBM/coeffs_poly.jl")
include("src/CLBM/collision_sym.jl")
include("src/CLBM/carleman_transferA.jl")
include("src/CLBM/carleman_transferA_ngrid.jl")
include("src/CLBM/LBM_const_subs.jl")
include("src/LBM/lbm_cons.jl")
include("src/LBM/lbm_const_sym.jl")
include("src/LBM/forcing.jl")
include("src/LBM/f_initial.jl")
include("src/CLBM/timeMarching.jl")
include("src/CLBM/CLBM_collision_test.jl")

println("=== CLBM PERFORMANCE ANALYSIS ===")
println("Configuration: ngrid=$ngrid, use_sparse=$use_sparse")

# Setup phase timing
println("\n1. SETUP PHASE TIMING:")

setup_time = @elapsed begin
    print("  Setting up LBM constants... ")
    const_time = @elapsed begin
        w, e, w_val, e_val = lbm_const_sym()
        global w_value = w_val
        global e_value = e_val
    end
    println("$(round(const_time, digits=3))s")
    
    print("  Computing collision operators... ")
    collision_time = @elapsed begin
        f, omega, u, rho = collision(Q, D, w, e, rho0, lTaylor, lorder2)
    end
    println("$(round(collision_time, digits=3))s")
    
    print("  Computing F coefficients... ")
    coeff_time = @elapsed begin
        global F1_ngrid, F2_ngrid, F3_ngrid = get_coeff_LBM_Fi_ngrid(poly_order, Q, f, omega, tau_value, ngrid)
    end
    println("$(round(coeff_time, digits=3))s")
end

println("  Total setup time: $(round(setup_time, digits=3))s")

# Matrix construction timing
println("\n2. MATRIX CONSTRUCTION TIMING:")

if use_sparse
    println("  Testing SPARSE matrix construction:")
    
    matrix_time = @elapsed begin
        C_sparse, bt_sparse, F0_sparse = carleman_C_sparse(Q, truncation_order, poly_order, f, omega, tau_value, force_factor, w_value, e_value)
    end
    println("  Sparse matrix construction: $(round(matrix_time, digits=3))s")
    
    matrix_size = size(C_sparse)
    nnz_count = nnz(C_sparse)
    sparsity = (1 - nnz_count / prod(matrix_size)) * 100
    
    println("  Matrix size: $(matrix_size[1])×$(matrix_size[2])")
    println("  Non-zeros: $nnz_count")
    println("  Sparsity: $(round(sparsity, digits=1))%")
    
else
    println("  Testing DENSE matrix construction:")
    
    matrix_time = @elapsed begin
        C_dense, bt_dense, F0_dense = carleman_C(Q, truncation_order, poly_order, f, omega, tau_value, force_factor, w_value, e_value)
    end
    println("  Dense matrix construction: $(round(matrix_time, digits=3))s")
    
    matrix_size = size(C_dense)
    println("  Matrix size: $(matrix_size[1])×$(matrix_size[2])")
    
    # Set the matrix for time marching test
    C_test = C_dense
end

# Time marching timing
println("\n3. TIME MARCHING TIMING:")

# Test with different numbers of time steps to see scaling
test_steps = [10, 50, 100]

for steps in test_steps
    print("  Testing $steps time steps... ")
    
    if use_sparse
        march_time = @elapsed begin
            fT, VT_f, VT = CLBM_collision_test_sparse(Q, omega, f, truncation_order, dt, tau_value, e_value, steps, false)
        end
    else
        march_time = @elapsed begin
            fT, VT_f, VT = CLBM_collision_test(Q, omega, f, C_test, truncation_order, dt, tau_value, e_value, steps, false)
        end
    end
    
    time_per_step = march_time / steps
    println("$(round(march_time, digits=3))s ($(round(time_per_step*1000, digits=1)) ms/step)")
end

# Detailed breakdown of time marching components
println("\n4. TIME MARCHING COMPONENT BREAKDOWN:")

steps = 50  # Use moderate number for detailed timing

if use_sparse
    println("  SPARSE time marching breakdown:")
    
    # Initial setup
    print("    Initial vector construction... ")
    init_time = @elapsed begin
        u0 = 0.1
        f_ini = f_ini_test(u0)
        V0 = carleman_V(f_ini, truncation_order)
        V0 = Float64.(V0)
    end
    println("$(round(init_time*1000, digits=1)) ms")
    
    # Matrix construction (already done above)
    print("    Sparse matrix construction... ")
    println("$(round(matrix_time*1000, digits=1)) ms (already measured)")
    
    # Time stepping loop
    print("    Time stepping loop ($steps steps)... ")
    
    # Recreate the setup for timing the loop only
    C_sparse, bt, F0 = carleman_C_sparse(Q, truncation_order, poly_order, f, omega, tau_value, force_factor, w_value, e_value)
    VT = zeros(length(V0), steps)
    VT[:, 1] = V0
    VT_f = zeros(Q, steps)
    VT_f[:, 1] = VT[1:Q, 1]
    uT = zeros(steps)
    _, uT[1] = lbm_u(e_value, VT_f[:, 1])
    
    omega_sub = LBM_const_subs(omega, tau_value)
    LB = lambdify(omega_sub * dt .+ f, f)
    fT = zeros(Q, steps)
    fT[:, 1] = f_ini
    
    loop_time = @elapsed begin
        for nt = 2:steps
            VT[:, nt] = (C_sparse * VT[:, nt - 1] + bt) .* dt .+ VT[:, nt - 1]
            _, uT[nt] = lbm_u(e_value, VT[1:Q, nt])
            fT_temp = LB(fT[1, nt-1], fT[2, nt-1], fT[3, nt-1]) + F0 * dt
            fT[:, nt] = fT_temp
        end
    end
    
    println("$(round(loop_time*1000, digits=1)) ms ($(round(loop_time*1000/steps, digits=2)) ms/step)")
    
    # Breakdown of single time step
    print("    Single matrix-vector multiply... ")
    test_v = VT[:, 1]
    matvec_time = @elapsed begin
        for i = 1:100  # Average over multiple ops
            result = C_sparse * test_v
        end
    end
    matvec_time /= 100
    println("$(round(matvec_time*1000, digits=3)) ms")
    
else
    println("  DENSE time marching breakdown:")
    
    # Similar breakdown for dense version
    print("    Initial vector construction... ")
    init_time = @elapsed begin
        u0 = 0.1
        f_ini = f_ini_test(u0)
        V0 = carleman_V(f_ini, truncation_order)
        V0 = Float64.(V0)
    end
    println("$(round(init_time*1000, digits=1)) ms")
    
    print("    Dense matrix construction... ")
    println("$(round(matrix_time*1000, digits=1)) ms (already measured)")
    
    # Similar time stepping analysis for dense
    print("    Time stepping loop ($steps steps)... ")
    
    C_dense = C_test
    _, bt, F0 = carleman_C(Q, truncation_order, poly_order, f, omega, tau_value, force_factor, w_value, e_value)
    VT = zeros(size(C_dense)[1], steps)
    VT[:, 1] = V0
    VT_f = zeros(Q, steps)
    VT_f[:, 1] = VT[1:Q, 1]
    uT = zeros(steps)
    _, uT[1] = lbm_u(e_value, VT_f[:, 1])
    
    omega_sub = LBM_const_subs(omega, tau_value)
    LB = lambdify(omega_sub * dt .+ f, f)
    fT = zeros(Q, steps)
    fT[:, 1] = f_ini
    
    loop_time = @elapsed begin
        for nt = 2:steps
            VT[:, nt] = (C_dense * VT[:, nt - 1] + bt) .* dt .+ VT[:, nt - 1]
            _, uT[nt] = lbm_u(e_value, VT[1:Q, nt])
            fT_temp = LB(fT[1, nt-1], fT[2, nt-1], fT[3, nt-1]) + F0 * dt
            fT[:, nt] = fT_temp
        end
    end
    
    println("$(round(loop_time*1000, digits=1)) ms ($(round(loop_time*1000/steps, digits=2)) ms/step)")
    
    print("    Single matrix-vector multiply... ")
    test_v = VT[:, 1]
    matvec_time = @elapsed begin
        for i = 1:100
            result = C_dense * test_v
        end
    end
    matvec_time /= 100
    println("$(round(matvec_time*1000, digits=3)) ms")
end

println("\n=== PERFORMANCE SUMMARY ===")
println("Setup time: $(round(setup_time, digits=3))s")
println("Matrix construction: $(round(matrix_time, digits=3))s")
total_time = setup_time + matrix_time
println("Total initialization: $(round(total_time, digits=3))s")

println("\nBottlenecks identified:")
if setup_time > matrix_time
    println("  1. Setup phase ($(round(setup_time/total_time*100, digits=1))% of init time)")
    println("  2. Matrix construction ($(round(matrix_time/total_time*100, digits=1))% of init time)")
else
    println("  1. Matrix construction ($(round(matrix_time/total_time*100, digits=1))% of init time)")
    println("  2. Setup phase ($(round(setup_time/total_time*100, digits=1))% of init time)")
end

println("\nFor repeated time stepping, each step takes ~$(round(loop_time*1000/steps, digits=2)) ms")
