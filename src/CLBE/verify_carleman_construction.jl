using Test
using LinearAlgebra
using SparseArrays
using Random

# Paper-aligned staged verification of the Carleman construction.
#
# This script is intended to verify the operator construction steps in
# Li_et_al.tex, especially:
#   - Eq. (lbm-carleman)
#   - Eq. (defiBij)
#   - Eq. (F_i_alpha)
#   - Eq. (diagonality-F)
#   - Eq. (npt-coll)
#
# The goal is to validate the algebraic construction of the Carleman matrix
# without relying on an expensive full D2Q9 CLBM time evolution.

include("clbe_config.jl")

l_sympy = true
QCFD_SRC = ENV["QCFD_SRC"]
QCFD_HOME = ENV["QCFD_HOME"]

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
include(QCFD_SRC * "LBM/lbm_const_sym.jl")
include(QCFD_SRC * "CLBE/timeMarching.jl")

function basis_vector(n, idx)
    v = zeros(Float64, n)
    v[idx] = 1.0
    return v
end

function state_kron(v, k)
    return Float64.(Kron_kth(v, k))
end

function point_velocity_index(alpha, vel, Q)
    return (alpha - 1) * Q + vel
end

function build_positive_probe_states(w_num; nstates=4, seed=1234)
    rng = MersenneTwister(seed)
    states = Vector{Vector{Float64}}()
    for _ = 1:nstates
        raw = abs.(Float64.(w_num) .+ 1.0e-2 .* randn(rng, length(w_num)))
        state = raw ./ sum(raw)
        push!(states, state)
    end
    return states
end

function d2q9_setup(; local_ngrid=1, local_rho0=1.0001, local_lTaylor=true)
    global Q = 9
    global D = 2
    global ngrid = local_ngrid
    global use_sparse = true
    global poly_order = 3
    global truncation_order = 3
    global rho0 = local_rho0
    global lTaylor = local_lTaylor
    global force_factor = 0.0

    w_sym, e_sym, w_num, e_num, a_sym, b_sym, c_sym, d_sym, a_val, b_val, c_val, d_val = lbm_const_sym(Q_local=Q, D_local=D)
    global w_value = Float64.(w_num)
    global e_value = e_num

    F1_single, F2_single, F3_single = numerical_carleman_coefficients(
        poly_order,
        Q,
        tau_value;
        w_value_input=w_value,
        e_value_input=e_value,
        rho_value_input=rho0,
        lTaylor_input=lTaylor,
        D_input=D,
    )

    return (
        w_sym=w_sym,
        e_sym=e_sym,
        w_num=w_num,
        e_num=e_num,
        a_sym=a_sym,
        b_sym=b_sym,
        c_sym=c_sym,
        d_sym=d_sym,
        a_val=a_val,
        b_val=b_val,
        c_val=c_val,
        d_val=d_val,
        F1_single=F1_single,
        F2_single=F2_single,
        F3_single=F3_single,
    )
end

old_settings = (
    Q=Q,
    D=D,
    ngrid=ngrid,
    use_sparse=use_sparse,
    poly_order=poly_order,
    truncation_order=truncation_order,
    rho0=rho0,
    lTaylor=lTaylor,
    force_factor=force_factor,
    w_value=copy(w_value),
    e_value=copy(e_value),
)

try
    setup = d2q9_setup(local_ngrid=1, local_rho0=1.0001, local_lTaylor=true)
    probe_states = build_positive_probe_states(setup.w_num)

    @testset "Paper-aligned Carleman construction checks (D2Q9)" begin
        @testset "Stage 1: Kronecker ordering for f^[k]" begin
            v = [1.0, 2.0, 3.0]

            expected_k2 = kron(v, v)
            expected_k3 = kron(v, kron(v, v))

            @test state_kron(v, 1) == v
            @test state_kron(v, 2) == expected_k2
            @test state_kron(v, 3) == expected_k3

            println("✅ Stage 1 passed: Kronecker ordering matches the f^[k] definition")
        end

        @testset "Stage 2: Eq. (lbm-carleman) polynomial collision decomposition" begin
            F1 = setup.F1_single
            F2 = setup.F2_single
            F3 = setup.F3_single

            for state in probe_states
                direct_rhs = numerical_collision_rhs(
                    state,
                    tau_value;
                    w_value_input=w_value,
                    e_value_input=e_value,
                    rho_value_input=rho0,
                    lTaylor_input=lTaylor,
                    D_input=D,
                )

                carleman_rhs = F1 * state + F2 * state_kron(state, 2) + F3 * state_kron(state, 3)

                @test carleman_rhs ≈ direct_rhs atol=1e-11 rtol=1e-11
            end

            println("✅ Stage 2 passed: F^(1), F^(2), F^(3) reproduce the D2Q9 polynomial collision RHS")
        end

        @testset "Stage 3: Eq. (defiBij) single-point transfer blocks" begin
            F1 = setup.F1_single
            F2 = setup.F2_single
            F3 = setup.F3_single

            A_2_2 = sum_Kron_kth_identity(F1, 2, Q)  # A_2^2 from F^(1)
            A_3_2 = sum_Kron_kth_identity(F2, 2, Q)  # A_3^2 from F^(2)
            A_4_2 = sum_Kron_kth_identity(F3, 2, Q)  # A_4^2 from F^(3)
            A_3_3 = sum_Kron_kth_identity(F1, 3, Q)  # A_3^3 from F^(1)

            for state in probe_states
                rhs1 = F1 * state
                rhs2 = F2 * state_kron(state, 2)
                rhs3 = F3 * state_kron(state, 3)

                expected_A_2_2 = kron(rhs1, state) + kron(state, rhs1)
                expected_A_3_2 = kron(rhs2, state) + kron(state, rhs2)
                expected_A_4_2 = kron(rhs3, state) + kron(state, rhs3)

                expected_A_3_3 =
                    kron(rhs1, state_kron(state, 2)) +
                    kron(state, kron(rhs1, state)) +
                    kron(state, kron(state, rhs1))

                @test A_2_2 * state_kron(state, 2) ≈ expected_A_2_2 atol=1e-11 rtol=1e-11
                @test A_3_2 * state_kron(state, 3) ≈ expected_A_3_2 atol=1e-11 rtol=1e-11
                @test A_4_2 * state_kron(state, 4) ≈ expected_A_4_2 atol=1e-11 rtol=1e-11
                @test A_3_3 * state_kron(state, 3) ≈ expected_A_3_3 atol=1e-11 rtol=1e-11
            end

            println("✅ Stage 3 passed: transfer blocks satisfy the Leibniz-lift construction")
        end

        @testset "Stage 4: Eq. (F_i_alpha) and Eq. (diagonality-F) n-point locality" begin
            global ngrid = 2

            F1_ngrid, F2_ngrid, F3_ngrid = get_coeff_LBM_Fi_ngrid(
                poly_order,
                Q,
                nothing,
                nothing,
                tau_value,
                ngrid;
                method=:numerical,
                w_value_input=w_value,
                e_value_input=e_value,
                rho_value_input=rho0,
                lTaylor_input=lTaylor,
                D_input=D,
            )

            F1 = setup.F1_single
            F2 = setup.F2_single
            F3 = setup.F3_single

            nphi = ngrid * Q

            idx_x1_v2 = point_velocity_index(1, 2, Q)
            idx_x1_v4 = point_velocity_index(1, 4, Q)
            idx_x2_v3 = point_velocity_index(2, 3, Q)
            idx_x2_v5 = point_velocity_index(2, 5, Q)

            e_x1_v2 = basis_vector(nphi, idx_x1_v2)
            e_x1_v4 = basis_vector(nphi, idx_x1_v4)
            e_x2_v3 = basis_vector(nphi, idx_x2_v3)
            e_x2_v5 = basis_vector(nphi, idx_x2_v5)

            out_F1_x1 = F1_ngrid * e_x1_v2
            @test out_F1_x1[1:Q] ≈ F1[:, 2] atol=1e-12 rtol=1e-12
            @test out_F1_x1[Q+1:2Q] ≈ zeros(Q) atol=1e-12 rtol=1e-12

            same_point_pair = kron(e_x1_v2, e_x1_v4)
            mixed_point_pair = kron(e_x1_v2, e_x2_v3)

            out_same_pair = F2_ngrid * same_point_pair
            out_mixed_pair = F2_ngrid * mixed_point_pair

            @test out_same_pair[1:Q] ≈ F2[:, index2(2, 4, Q)] atol=1e-12 rtol=1e-12
            @test out_same_pair[Q+1:2Q] ≈ zeros(Q) atol=1e-12 rtol=1e-12
            @test out_mixed_pair ≈ zeros(nphi) atol=1e-12 rtol=1e-12

            same_point_triple = kron(e_x2_v3, kron(e_x2_v5, e_x2_v3))
            mixed_point_triple = kron(e_x1_v2, kron(e_x1_v4, e_x2_v5))

            out_same_triple = F3_ngrid * same_point_triple
            out_mixed_triple = F3_ngrid * mixed_point_triple

            @test out_same_triple[1:Q] ≈ zeros(Q) atol=1e-12 rtol=1e-12
            @test out_same_triple[Q+1:2Q] ≈ F3[:, index3(3, 5, 3, Q)] atol=1e-12 rtol=1e-12
            @test out_mixed_triple ≈ zeros(nphi) atol=1e-12 rtol=1e-12

            println("✅ Stage 4 passed: n-point lifted collision operators preserve spatial locality")
        end

        @testset "Stage 5: Eq. (npt-coll) assembled n-point collision Carleman matrix" begin
            global ngrid = 2

            global F1_ngrid, F2_ngrid, F3_ngrid = get_coeff_LBM_Fi_ngrid(
                poly_order,
                Q,
                nothing,
                nothing,
                tau_value,
                ngrid;
                method=:numerical,
                w_value_input=w_value,
                e_value_input=e_value,
                rho_value_input=rho0,
                lTaylor_input=lTaylor,
                D_input=D,
            )

            phi_probe = vcat(probe_states[1], probe_states[2])
            V0 = Float64.(carleman_V(phi_probe, truncation_order))

            C_sparse, bt_sparse, _ = carleman_C_sparse(
                Q,
                truncation_order,
                poly_order,
                nothing,
                nothing,
                tau_value,
                force_factor,
                w_value,
                e_value,
            )

            direct_collision_rhs =
                F1_ngrid * phi_probe +
                F2_ngrid * state_kron(phi_probe, 2) +
                F3_ngrid * state_kron(phi_probe, 3)

            carleman_collision_rhs = Array(C_sparse * V0 + bt_sparse)

            @test size(C_sparse, 1) == carleman_C_dim(Q, truncation_order, ngrid)
            @test size(C_sparse, 2) == carleman_C_dim(Q, truncation_order, ngrid)
            @test carleman_collision_rhs[1:ngrid*Q] ≈ direct_collision_rhs atol=1e-11 rtol=1e-11

            println("✅ Stage 5 passed: assembled n-point collision Carleman matrix reproduces the direct collision RHS")
        end
    end
finally
    global Q = old_settings.Q
    global D = old_settings.D
    global ngrid = old_settings.ngrid
    global use_sparse = old_settings.use_sparse
    global poly_order = old_settings.poly_order
    global truncation_order = old_settings.truncation_order
    global rho0 = old_settings.rho0
    global lTaylor = old_settings.lTaylor
    global force_factor = old_settings.force_factor
    global w_value = old_settings.w_value
    global e_value = old_settings.e_value
end

println("\n🎉 D2Q9 paper-aligned Carleman construction checks completed successfully!")