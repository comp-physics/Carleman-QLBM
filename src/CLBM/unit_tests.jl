# Unit tests for individual CLBM functions
using Test
using SparseArrays
using LinearAlgebra

# Load configuration and functions
include("clbm_config.jl")

# Set up symbolic computation
l_sympy = true
QCFD_SRC = ENV["QCFD_SRC"]
QCFD_HOME = ENV["QCFD_HOME"]

# Include necessary modules
if l_sympy
    using SymPy
    include(QCFD_SRC * "CLBM/coeffs_poly.jl")
else
    using Symbolics
end

include(QCFD_SRC * "CLBM/collision_sym.jl")
include(QCFD_SRC * "CLBM/carleman_transferA.jl")
include(QCFD_SRC * "CLBM/carleman_transferA_ngrid.jl") 
include(QCFD_SRC * "CLBM/LBM_const_subs.jl")
include(QCFD_SRC * "LBM/lbm_cons.jl")
include(QCFD_SRC * "LBM/lbm_const_sym.jl")
include(QCFD_SRC * "LBM/forcing.jl")
include(QCFD_SRC * "LBM/f_initial.jl")
include(QCFD_SRC * "CLBM/timeMarching.jl")

@testset "CLBM Unit Tests" begin
    
    @testset "Configuration Tests" begin
        @test Q == 3
        @test D == 1
        @test truncation_order == 3
        @test poly_order == 3
        @test isa(tau_value, Float64)
        @test n_time > 0
    end
    
    @testset "Sparse Kronecker Functions" begin
        # Test small matrices
        A = [1.0 2.0; 3.0 4.0]
        
        # Test Kron_kth_sparse
        A_sparse_1 = Kron_kth_sparse(A, 1)
        @test A_sparse_1 ≈ sparse(A)
        
        A_sparse_2 = Kron_kth_sparse(A, 2)
        @test A_sparse_2 ≈ sparse(kron(A, A))
        
        # Test with identity matrices
        I3 = Matrix(1.0I, 3, 3)
        I3_sparse = Kron_kth_sparse(I3, 2)
        @test I3_sparse ≈ sparse(kron(I3, I3))
        
        println("✅ Sparse Kronecker functions work correctly")
    end
    
    @testset "Matrix Dimensions" begin
        # Test Carleman matrix dimension calculation
        C_dim = carleman_C_dim(Q, truncation_order, ngrid)
        @test C_dim > 0
        @test isa(C_dim, Int)
        
        # For Q=3, truncation_order=3, ngrid=1, expect specific dimension
        if Q == 3 && truncation_order == 3 && ngrid == 1
            @test C_dim == 39  # 3 + 9 + 27 = 39
        end
        
        println("✅ Matrix dimensions calculated correctly")
    end
    
    @testset "LBM Setup" begin
        # Test LBM constants
        w, e, w_val, e_val = lbm_const_sym()
        @test length(w_val) == Q
        @test length(e_val) == Q
        
        # Test collision operator generation
        f, omega, u, rho = collision(Q, D, w, e, rho0, lTaylor, lorder2)
        @test length(f) == Q
        @test length(omega) == Q
        
        println("✅ LBM setup functions work correctly")
    end
    
    @testset "Sparse vs Dense Matrix Construction" begin
        # Set up LBM
        w, e, w_val, e_val = lbm_const_sym()
        global w_value = w_val
        global e_value = e_val
        f, omega, u, rho = collision(Q, D, w, e, rho0, lTaylor, lorder2)
        
        # Initialize F matrices
        global F1_ngrid, F2_ngrid, F3_ngrid = get_coeff_LBM_Fi_ngrid(poly_order, Q, f, omega, tau_value, ngrid)
        
        # Test that we can build both matrices
        C_dense, bt_dense, F0_dense = carleman_C(Q, truncation_order, poly_order, f, omega, tau_value, force_factor, w_val, e_val)
        C_sparse, bt_sparse, F0_sparse = carleman_C_sparse(Q, truncation_order, poly_order, f, omega, tau_value, force_factor, w_val, e_val)
        
        # Basic checks
        @test size(C_dense) == size(C_sparse)
        @test length(bt_dense) == length(bt_sparse)
        @test F0_dense ≈ F0_sparse
        
        # Check sparsity
        sparsity = (length(C_sparse.nzval) / (size(C_sparse, 1) * size(C_sparse, 2))) * 100
        @test sparsity < 50  # Should be quite sparse
        
        println("✅ Sparse and dense matrix construction both work")
    end
    
    @testset "Initial Conditions" begin
        # Test initial condition generation
        u0 = 0.1
        f_ini = f_ini_test(u0)
        @test length(f_ini) == Q
        @test all(f_ini .>= 0)  # Distribution functions should be non-negative
        
        # Test Carleman vector generation
        V = carleman_V(f_ini, truncation_order)
        expected_length = carleman_C_dim(Q, truncation_order, ngrid)
        @test length(V) == expected_length
        
        println("✅ Initial conditions generated correctly")
    end
end

println("\n🎉 All unit tests completed!")
