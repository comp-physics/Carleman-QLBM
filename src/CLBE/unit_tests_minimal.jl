# Minimal unit tests for CI environments (no plotting dependencies)
using Test
using SparseArrays
using LinearAlgebra

# Load configuration and functions
include("clbe_config.jl")

# Set up symbolic computation
l_sympy = true
QCFD_SRC = ENV["QCFD_SRC"]
QCFD_HOME = ENV["QCFD_HOME"]

# Include necessary modules (skip plotting)
if l_sympy
    using SymPy
    include(QCFD_SRC * "CLBE/coeffs_poly.jl")
else
    using Symbolics
end

include(QCFD_SRC * "CLBE/collision_sym.jl")
include(QCFD_SRC * "CLBE/carleman_transferA.jl")
include(QCFD_SRC * "CLBE/carleman_transferA_ngrid.jl") 
include(QCFD_SRC * "CLBE/LBM_const_subs.jl")
include(QCFD_SRC * "LBM/lbm_cons.jl")
include(QCFD_SRC * "LBM/lbm_const_sym.jl")
include(QCFD_SRC * "LBM/forcing.jl")
include(QCFD_SRC * "LBM/f_initial.jl")
include(QCFD_SRC * "CLBE/streaming_Carleman.jl")
include(QCFD_SRC * "CLBE/timeMarching.jl")
include(QCFD_SRC * "LBE/direct_LBE.jl")

@testset "CLBE Minimal Tests" begin
    
    @testset "Configuration Tests" begin
        @test Q == 3
        @test D == 1
        @test truncation_order == 3
        @test poly_order == 3
        @test isa(tau_value, Float64)
        @test n_time > 0
        println("✅ Configuration validated")
    end
    
    @testset "Basic Function Loading" begin
        # Test that key functions can be called without errors
        w, e, w_val, e_val = lbm_const_sym()
        @test length(w_val) == Q
        @test length(e_val) == Q
        
        f, omega, u, rho = collision(Q, D, w, e, rho0, lTaylor, lorder2)
        @test length(f) == Q
        @test length(omega) == Q
        
        println("✅ LBM functions load correctly")
    end
    
    @testset "Matrix Dimension Calculation" begin
        C_dim = carleman_C_dim(Q, truncation_order, ngrid)
        @test C_dim > 0
        @test isa(C_dim, Int)
        
        # For Q=3, truncation_order=3, ngrid=1, expect specific dimension
        if Q == 3 && truncation_order == 3 && ngrid == 1
            @test C_dim == 39  # 3 + 9 + 27 = 39
        end
        
        println("✅ Matrix dimensions calculated correctly")
    end
    
    @testset "Sparse Kronecker Functions" begin
        # Test small matrices to avoid memory issues
        A = [1.0 2.0; 3.0 4.0]
        
        # Test basic sparse conversion
        A_sparse = sparse(A)
        @test A_sparse ≈ A
        
        # Test Kronecker product with itself
        A_kron = kron(A, A)
        A_kron_sparse = sparse(A_kron)
        @test A_kron_sparse ≈ A_kron
        
        println("✅ Sparse operations work correctly")
    end
    
    @testset "Initial Conditions" begin
        # Test initial condition generation
        u0 = 0.1
        f_ini = f_ini_test(u0)
        @test length(f_ini) == Q
        @test all(f_ini .>= 0)  # Distribution functions should be non-negative

        sinusoidal_velocity_profile = d1q3_sinusoidal_velocity_profile(6; u_ini=u0)
        @test length(sinusoidal_velocity_profile) == 6
        @test isapprox(sinusoidal_velocity_profile[1], u0 * sin(2 * pi / 6), atol=1e-12)

        phi_sinusoidal = d1q3_multigrid_initial_condition(6; initial_condition=:sinusoidal, u_ini=u0)
        @test length(phi_sinusoidal) == Q * 6
        @test all(phi_sinusoidal .>= 0)
        
        # Test Carleman vector generation
        V = carleman_V(f_ini, truncation_order)
        expected_length = carleman_C_dim(Q, truncation_order, ngrid)
        @test length(V) == expected_length
        
        println("✅ Initial conditions work correctly")
    end

    @testset "Numerical coefficient generation matches symbolic D1Q3" begin
        old_ngrid = ngrid

        try
            global ngrid = 3

            w, e, w_val, e_val = lbm_const_sym()
            global w_value = w_val
            global e_value = e_val

            f, omega, u, rho = collision(Q, D, w, e, rho0, lTaylor, lorder2)

            F1_symbolic, F2_symbolic, F3_symbolic = get_coeff_LBM_Fi_ngrid(
                poly_order,
                Q,
                f,
                omega,
                tau_value,
                ngrid;
                method=:symbolic,
            )

            F1_numeric, F2_numeric, F3_numeric = get_coeff_LBM_Fi_ngrid(
                poly_order,
                Q,
                f,
                omega,
                tau_value,
                ngrid;
                method=:numerical,
                w_value_input=w_val,
                e_value_input=e_val,
                rho_value_input=rho0,
                lTaylor_input=lTaylor,
                D_input=D,
            )

            @test F1_numeric ≈ F1_symbolic atol=1e-12 rtol=1e-12
            @test F2_numeric ≈ F2_symbolic atol=1e-12 rtol=1e-12
            @test F3_numeric ≈ F3_symbolic atol=1e-12 rtol=1e-12

            phi_ini = vcat(
                f_ini_test(0.12),
                f_ini_test(0.00),
                f_ini_test(-0.08),
            )
            S_lbm, _ = streaming_operator_D1Q3_interleaved(ngrid, 1)

            rhs_symbolic = direct_lbe_rhs_ngrid(phi_ini, S_lbm, F1_symbolic, F2_symbolic, F3_symbolic)
            rhs_numeric = direct_lbe_rhs_ngrid(phi_ini, S_lbm, F1_numeric, F2_numeric, F3_numeric)
            @test rhs_numeric ≈ rhs_symbolic atol=1e-12 rtol=1e-12

            global F1_ngrid, F2_ngrid, F3_ngrid = F1_numeric, F2_numeric, F3_numeric
            phiT_numeric, VT_numeric = timeMarching_state_CLBM_sparse(omega, f, tau_value, Q, truncation_order, dt, phi_ini, 6; S_lbm=S_lbm)

            global F1_ngrid, F2_ngrid, F3_ngrid = F1_symbolic, F2_symbolic, F3_symbolic
            phiT_symbolic, VT_symbolic = timeMarching_state_CLBM_sparse(omega, f, tau_value, Q, truncation_order, dt, phi_ini, 6; S_lbm=S_lbm)

            @test phiT_numeric ≈ phiT_symbolic atol=1e-12 rtol=1e-12
            @test VT_numeric ≈ VT_symbolic atol=1e-12 rtol=1e-12

            println("✅ Numerical D1Q3 coefficients reproduce symbolic coefficients and dynamics")
        finally
            global ngrid = old_ngrid
        end
    end

    @testset "ngrid=2 Regression" begin
        old_ngrid = ngrid
        old_use_sparse = use_sparse
        old_n_time = n_time

        try
            global ngrid = 2
            global use_sparse = true
            local_n_time = 10

            w, e, w_val, e_val = lbm_const_sym()
            global w_value = w_val
            global e_value = e_val

            f, omega, u, rho = collision(Q, D, w, e, rho0, lTaylor, lorder2)
            global F1_ngrid, F2_ngrid, F3_ngrid = get_coeff_LBM_Fi_ngrid(poly_order, Q, f, omega, tau_value, ngrid)

            f_ini = f_ini_test(u0)
            phi_ini = get_phi(f_ini, ngrid)

            F1_single = F_carlemanOrder_collision(Q, 1, f, omega, tau_value)
            F2_single = F_carlemanOrder_collision(Q, 2, f, omega, tau_value)
            F3_single = F_carlemanOrder_collision(Q, 3, f, omega, tau_value)
            single_rhs = F1_single * f_ini + F2_single * kron(f_ini, f_ini) + F3_single * kron(f_ini, kron(f_ini, f_ini))
            multi_rhs = F1_ngrid * phi_ini + F2_ngrid * kron(phi_ini, phi_ini) + F3_ngrid * kron(phi_ini, kron(phi_ini, phi_ini))

            @test multi_rhs[1:Q] ≈ single_rhs atol=1e-12 rtol=1e-12
            @test multi_rhs[Q+1:2Q] ≈ single_rhs atol=1e-12 rtol=1e-12

            VT_f, VT, uT, fT = timeMarching_collision_CLBM_sparse(
                omega, f, tau_value, Q, truncation_order, e_val, dt, f_ini, local_n_time, false
            )

            @test VT_f[:, end] ≈ fT[:, end] atol=1e-12 rtol=1e-12
            @test minimum(VT_f) ≥ -1e-12

            println("✅ ngrid=2 periodic CLBM matches LBM")
        finally
            global ngrid = old_ngrid
            global use_sparse = old_use_sparse
            global n_time = old_n_time
        end
    end

    @testset "ngrid=3 Nonuniform Streaming+Collision" begin
        old_ngrid = ngrid
        old_use_sparse = use_sparse

        try
            global ngrid = 3
            global use_sparse = true

            w, e, w_val, e_val = lbm_const_sym()
            global w_value = w_val
            global e_value = e_val

            f, omega, u, rho = collision(Q, D, w, e, rho0, lTaylor, lorder2)
            global F1_ngrid, F2_ngrid, F3_ngrid = get_coeff_LBM_Fi_ngrid(poly_order, Q, f, omega, tau_value, ngrid)

            phi_ini = vcat(
                f_ini_test(0.12),
                f_ini_test(0.00),
                f_ini_test(-0.08),
            )

            S_lbm, _ = streaming_operator_D1Q3_interleaved(ngrid, 1)
            direct_rhs = -S_lbm * phi_ini +
                F1_ngrid * phi_ini +
                F2_ngrid * kron(phi_ini, phi_ini) +
                F3_ngrid * kron(phi_ini, kron(phi_ini, phi_ini))

            C_sparse, bt_sparse, _ = carleman_C_sparse(Q, truncation_order, poly_order, f, omega, tau_value, force_factor, w_val, e_val)
            S_sparse = build_streaming_carleman_operator_sparse(Q, truncation_order, poly_order, ngrid)
            V0 = Float64.(carleman_V(phi_ini, truncation_order))
            carleman_rhs = Array((C_sparse - S_sparse) * V0 + bt_sparse)

            @test carleman_rhs[1:Q*ngrid] ≈ direct_rhs atol=1e-12 rtol=1e-12

            direct_next = phi_ini + dt * direct_rhs
            carleman_next = phi_ini + dt * carleman_rhs[1:Q*ngrid]

            @test carleman_next ≈ direct_next atol=1e-12 rtol=1e-12
            @test all(isfinite.(carleman_next))

            println("✅ ngrid=3 nonuniform centered-difference streaming + collision matches direct n-point LBE")
        finally
            global ngrid = old_ngrid
            global use_sparse = old_use_sparse
        end
    end
end

println("\n🎉 Minimal unit tests completed successfully!")
