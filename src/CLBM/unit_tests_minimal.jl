# Minimal unit tests for CI environments (no plotting dependencies)
using Test
using SparseArrays
using LinearAlgebra

# Load configuration and functions
include("clbm_config.jl")

# Set up symbolic computation
l_sympy = true
QCFD_SRC = ENV["QCFD_SRC"]
QCFD_HOME = ENV["QCFD_HOME"]

# Include necessary modules (skip plotting)
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

@testset "CLBM Minimal Tests" begin
    
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
        
        # Test Carleman vector generation
        V = carleman_V(f_ini, truncation_order)
        expected_length = carleman_C_dim(Q, truncation_order, ngrid)
        @test length(V) == expected_length
        
        println("✅ Initial conditions work correctly")
    end
end

println("\n🎉 Minimal unit tests completed successfully!")
