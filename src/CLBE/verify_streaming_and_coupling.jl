using Test
using LinearAlgebra
using SparseArrays
using Random

# Staged verification of the streaming operator and collision+streaming coupling.
#
# This script follows the construction logic from Li_et_al.tex, especially the
# "Carleman-linearized LBE on spatial discretization grids" subsection:
#   - streaming operator on the interleaved spatial state φ(x),
#   - lifted streaming blocks B_i^i,
#   - assembled streaming operator C_s^(k),
#   - assembled collision operator C_c^(k),
#   - coupled operator C^(k) = C_s^(k) + C_c^(k).
#
# The validation is done in two layers:
#   1. D1Q3 on a small spatial grid, where direct nonlinear spatial LBE and
#      truncated CLBM can be compared cheaply and explicitly.
#   2. D2Q9 on tiny periodic grids, where the raw periodic stencil, lifted
#      streaming operator, and first-level coupled RHS are verified.

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
include(QCFD_SRC * "LBM/f_initial.jl")
include(QCFD_SRC * "CLBE/streaming_Carleman.jl")
include(QCFD_SRC * "CLBE/timeMarching.jl")
include(QCFD_SRC * "LBE/direct_LBE.jl")

function basis_vector(n, idx)
    v = zeros(Float64, n)
    v[idx] = 1.0
    return v
end

function state_kron(v, k)
    return Float64.(Kron_kth(v, k))
end

function d1q3_interleaved_index(pos, vel)
    return (pos - 1) * 3 + vel
end

function d2q9_interleaved_index(i, j, vel, nx)
    return ((j - 1) * nx + (i - 1)) * 9 + vel
end

function build_sitewise_positive_state(reference_weights, nsites; seed=1234, perturbation=1.0e-2)
    rng = MersenneTwister(seed)
    phi = Float64[]
    for _ = 1:nsites
        raw = abs.(Float64.(reference_weights) .+ perturbation .* randn(rng, length(reference_weights)))
        append!(phi, raw ./ sum(raw))
    end
    return phi
end

function current_d1q3_streaming_velocities()
    # This reflects the current implementation in streaming_operator_D1Q3_interleaved:
    # [rest, +x, -x]
    return [0.0, 1.0, -1.0]
end

function expected_d1q3_streaming_row(nx, hx, pos, vel)
    row = zeros(Float64, 3 * nx)
    ex = current_d1q3_streaming_velocities()[vel]

    if ex != 0.0
        left_pos = pos == 1 ? nx : pos - 1
        right_pos = pos == nx ? 1 : pos + 1
        row[d1q3_interleaved_index(left_pos, vel)] += -ex / (2 * hx)
        row[d1q3_interleaved_index(right_pos, vel)] += ex / (2 * hx)
    end

    return row
end

function d2q9_periodic_velocities()
    return [
        [0.0, 0.0],
        [1.0, 0.0],
        [-1.0, 0.0],
        [0.0, 1.0],
        [0.0, -1.0],
        [1.0, 1.0],
        [1.0, -1.0],
        [-1.0, 1.0],
        [-1.0, -1.0],
    ]
end

function streaming_operator_D2Q9_interleaved_periodic_test(nx, ny, hx, hy)
    e = d2q9_periodic_velocities()
    n_velocities = 9
    n_total = n_velocities * nx * ny

    I_idx = Int[]
    J_idx = Int[]
    V_vals = Float64[]

    for j in 1:ny
        for i in 1:nx
            for vel in 1:n_velocities
                row_idx = d2q9_interleaved_index(i, j, vel, nx)
                ex, ey = e[vel]

                if ex == 0.0 && ey == 0.0
                    push!(I_idx, row_idx)
                    push!(J_idx, row_idx)
                    push!(V_vals, 0.0)
                    continue
                end

                if ex != 0.0
                    left_i = i == 1 ? nx : i - 1
                    right_i = i == nx ? 1 : i + 1

                    push!(I_idx, row_idx)
                    push!(J_idx, d2q9_interleaved_index(left_i, j, vel, nx))
                    push!(V_vals, -ex / (2 * hx))

                    push!(I_idx, row_idx)
                    push!(J_idx, d2q9_interleaved_index(right_i, j, vel, nx))
                    push!(V_vals, ex / (2 * hx))
                end

                if ey != 0.0
                    bottom_j = j == 1 ? ny : j - 1
                    top_j = j == ny ? 1 : j + 1

                    push!(I_idx, row_idx)
                    push!(J_idx, d2q9_interleaved_index(i, bottom_j, vel, nx))
                    push!(V_vals, -ey / (2 * hy))

                    push!(I_idx, row_idx)
                    push!(J_idx, d2q9_interleaved_index(i, top_j, vel, nx))
                    push!(V_vals, ey / (2 * hy))
                end
            end
        end
    end

    return sparse(I_idx, J_idx, V_vals, n_total, n_total), e
end

function expected_d2q9_periodic_streaming_row(nx, ny, hx, hy, i, j, vel)
    row = zeros(Float64, 9 * nx * ny)
    ex, ey = d2q9_periodic_velocities()[vel]

    if ex != 0.0
        left_i = i == 1 ? nx : i - 1
        right_i = i == nx ? 1 : i + 1
        row[d2q9_interleaved_index(left_i, j, vel, nx)] += -ex / (2 * hx)
        row[d2q9_interleaved_index(right_i, j, vel, nx)] += ex / (2 * hx)
    end

    if ey != 0.0
        bottom_j = j == 1 ? ny : j - 1
        top_j = j == ny ? 1 : j + 1
        row[d2q9_interleaved_index(i, bottom_j, vel, nx)] += -ey / (2 * hy)
        row[d2q9_interleaved_index(i, top_j, vel, nx)] += ey / (2 * hy)
    end

    return row
end

function d1q3_setup(; local_ngrid=3, local_rho0=1.0001, local_lTaylor=lTaylor)
    global Q = 3
    global D = 1
    global ngrid = local_ngrid
    global use_sparse = true
    global poly_order = 3
    global truncation_order = 3
    global rho0 = local_rho0
    global lTaylor = local_lTaylor
    global force_factor = 0.0

    w_sym, e_sym, w_num, e_num, a_sym, b_sym, c_sym, d_sym, a_val, b_val, c_val, d_val = lbm_const_sym(Q_local=Q, D_local=D)
    global w_value = Float64.(w_num)
    global e_value = Float64.(e_num)

    f, omega, u, rho = collision(Q, D, w_sym, e_sym, rho0, lTaylor, lorder2)

    F1_local, F2_local, F3_local = get_coeff_LBM_Fi_ngrid(
        poly_order,
        Q,
        f,
        omega,
        tau_value,
        ngrid;
        method=:numerical,
        w_value_input=w_value,
        e_value_input=e_value,
        rho_value_input=rho0,
        lTaylor_input=lTaylor,
        D_input=D,
    )

    global F1_ngrid = F1_local
    global F2_ngrid = F2_local
    global F3_ngrid = F3_local

    S_lbm, _ = streaming_operator_D1Q3_interleaved(ngrid, 1.0)

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
        f=f,
        omega=omega,
        S_lbm=S_lbm,
        F1_ngrid=F1_local,
        F2_ngrid=F2_local,
        F3_ngrid=F3_local,
    )
end

function d2q9_setup(; local_nx=2, local_ny=2, local_rho0=1.0001, local_lTaylor=true)
    global Q = 9
    global D = 2
    global ngrid = local_nx * local_ny
    global use_sparse = true
    global poly_order = 3
    global truncation_order = 3
    global rho0 = local_rho0
    global lTaylor = local_lTaylor
    global force_factor = 0.0

    w_sym, e_sym, w_num, e_num, a_sym, b_sym, c_sym, d_sym, a_val, b_val, c_val, d_val = lbm_const_sym(Q_local=Q, D_local=D)
    global w_value = Float64.(w_num)
    global e_value = e_num

    F1_local, F2_local, F3_local = get_coeff_LBM_Fi_ngrid(
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

    global F1_ngrid = F1_local
    global F2_ngrid = F2_local
    global F3_ngrid = F3_local

    S_lbm, _ = streaming_operator_D2Q9_interleaved_periodic_test(local_nx, local_ny, 1.0, 1.0)

    return (
        nx=local_nx,
        ny=local_ny,
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
        S_lbm=S_lbm,
        F1_ngrid=F1_local,
        F2_ngrid=F2_local,
        F3_ngrid=F3_local,
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
    F1_ngrid=@isdefined(F1_ngrid) ? F1_ngrid : nothing,
    F2_ngrid=@isdefined(F2_ngrid) ? F2_ngrid : nothing,
    F3_ngrid=@isdefined(F3_ngrid) ? F3_ngrid : nothing,
)

try
    @testset "Streaming and coupled Carleman validation" begin
        @testset "D1Q3 staged streaming + coupling checks" begin
            setup = d1q3_setup(local_ngrid=3, local_rho0=1.0001, local_lTaylor=lTaylor)
            nphi = Q * ngrid

            phi_probe = vcat(
                f_ini_test(0.12),
                f_ini_test(0.00),
                f_ini_test(-0.08),
            )

            @testset "Stage 1: interleaved spatial state layout" begin
                @test phi_probe[1:3] == f_ini_test(0.12)
                @test phi_probe[4:6] == f_ini_test(0.00)
                @test phi_probe[7:9] == f_ini_test(-0.08)

                println("✅ D1Q3 Stage 1 passed: current interleaved spatial state layout is consistent site-by-site")
            end

            @testset "Stage 2: base centered-difference streaming stencil" begin
                for pos in 1:ngrid
                    for vel in 1:Q
                        row_idx = d1q3_interleaved_index(pos, vel)
                        actual_row = vec(Array(setup.S_lbm[row_idx, :]))
                        expected_row = expected_d1q3_streaming_row(ngrid, 1.0, pos, vel)
                        @test actual_row ≈ expected_row atol=1e-12 rtol=1e-12
                    end
                end

                println("✅ D1Q3 Stage 2 passed: raw interleaved streaming matrix matches the periodic centered-difference stencil")
            end

            @testset "Stage 3: lifted streaming operator satisfies the product rule" begin
                S_sparse = build_streaming_carleman_operator_sparse(Q, truncation_order, poly_order, ngrid; S_lbm=setup.S_lbm)
                block22, _ = carleman_C_block_dim(Q, 2, 2, 0)
                block33, _ = carleman_C_block_dim(Q, 3, 3, 0)

                Sphi = setup.S_lbm * phi_probe
                phi2 = state_kron(phi_probe, 2)
                phi3 = state_kron(phi_probe, 3)

                expected_order2 = kron(Sphi, phi_probe) + kron(phi_probe, Sphi)
                expected_order3 =
                    kron(Sphi, kron(phi_probe, phi_probe)) +
                    kron(phi_probe, kron(Sphi, phi_probe)) +
                    kron(phi_probe, kron(phi_probe, Sphi))

                actual_order2 = Array(S_sparse[block22, block22] * phi2)
                actual_order3 = Array(S_sparse[block33, block33] * phi3)

                @test actual_order2 ≈ expected_order2 atol=1e-12 rtol=1e-12
                @test actual_order3 ≈ expected_order3 atol=1e-12 rtol=1e-12

                println("✅ D1Q3 Stage 3 passed: lifted streaming blocks reproduce the product-rule action on φ^[2] and φ^[3]")
            end

            @testset "Stage 4: streaming-only first level" begin
                V0 = Float64.(carleman_V(phi_probe, truncation_order))
                S_sparse = build_streaming_carleman_operator_sparse(Q, truncation_order, poly_order, ngrid; S_lbm=setup.S_lbm)

                direct_streaming_rhs = -setup.S_lbm * phi_probe
                lifted_streaming_rhs = Array(-S_sparse * V0)

                @test lifted_streaming_rhs[1:nphi] ≈ direct_streaming_rhs atol=1e-12 rtol=1e-12

                println("✅ D1Q3 Stage 4 passed: first Carleman level reproduces streaming-only RHS")
            end

            @testset "Stage 5: coupled collision + streaming RHS" begin
                V0 = Float64.(carleman_V(phi_probe, truncation_order))
                C_sparse, bt_sparse, _ = carleman_C_sparse(Q, truncation_order, poly_order, setup.f, setup.omega, tau_value, force_factor, w_value, e_value)
                S_sparse = build_streaming_carleman_operator_sparse(Q, truncation_order, poly_order, ngrid; S_lbm=setup.S_lbm)

                direct_rhs = direct_lbe_rhs_ngrid(phi_probe, setup.S_lbm, setup.F1_ngrid, setup.F2_ngrid, setup.F3_ngrid)
                lifted_rhs = Array((C_sparse - S_sparse) * V0 + bt_sparse)

                @test lifted_rhs[1:nphi] ≈ direct_rhs atol=1e-12 rtol=1e-12

                println("✅ D1Q3 Stage 5 passed: assembled coupled Carleman operator reproduces the direct nonlinear spatial LBE RHS")
            end

            @testset "Stage 6: one-step explicit Euler update" begin
                V0 = Float64.(carleman_V(phi_probe, truncation_order))
                C_sparse, bt_sparse, _ = carleman_C_sparse(Q, truncation_order, poly_order, setup.f, setup.omega, tau_value, force_factor, w_value, e_value)
                S_sparse = build_streaming_carleman_operator_sparse(Q, truncation_order, poly_order, ngrid; S_lbm=setup.S_lbm)

                direct_rhs = direct_lbe_rhs_ngrid(phi_probe, setup.S_lbm, setup.F1_ngrid, setup.F2_ngrid, setup.F3_ngrid)
                lifted_rhs = Array((C_sparse - S_sparse) * V0 + bt_sparse)

                direct_next = phi_probe + dt * direct_rhs
                lifted_next = phi_probe + dt * lifted_rhs[1:nphi]

                @test lifted_next ≈ direct_next atol=1e-12 rtol=1e-12

                println("✅ D1Q3 Stage 6 passed: one explicit Euler step agrees between direct spatial LBE and CLBM first level")
            end

            @testset "Stage 7: short-time trajectory entry check" begin
                local_n_time = 6

                direct_history = timeMarching_direct_LBE_ngrid(phi_probe, dt, local_n_time, setup.F1_ngrid, setup.F2_ngrid, setup.F3_ngrid; S_lbm=setup.S_lbm)
                clbm_history, _ = timeMarching_state_CLBM_sparse(setup.omega, setup.f, tau_value, Q, truncation_order, dt, phi_probe, local_n_time; S_lbm=setup.S_lbm, nspatial=ngrid)

                @test clbm_history[:, 1] ≈ direct_history[:, 1] atol=1e-12 rtol=1e-12
                @test clbm_history[:, 2] ≈ direct_history[:, 2] atol=1e-12 rtol=1e-12
                @test all(isfinite.(clbm_history))
                @test all(isfinite.(direct_history))

                println("✅ D1Q3 Stage 7 passed: short-time direct and CLBM histories agree at initialization and first update, with finite trajectories")
            end
        end

        @testset "D2Q9 staged streaming + coupling checks" begin
            @testset "Stage 1: raw periodic D2Q9 streaming stencil" begin
                nx = 3
                ny = 3
                S_lbm, _ = streaming_operator_D2Q9_interleaved_periodic_test(nx, ny, 1.0, 1.0)

                for (i, j) in ((1, 1), (2, 2), (3, 3))
                    for vel in 1:9
                        row_idx = d2q9_interleaved_index(i, j, vel, nx)
                        actual_row = vec(Array(S_lbm[row_idx, :]))
                        expected_row = expected_d2q9_periodic_streaming_row(nx, ny, 1.0, 1.0, i, j, vel)
                        @test actual_row ≈ expected_row atol=1e-12 rtol=1e-12
                    end
                end

                println("✅ D2Q9 Stage 1 passed: periodic interleaved streaming matrix matches the expected x/y centered-difference stencil")
            end

            setup = d2q9_setup(local_nx=2, local_ny=2, local_rho0=1.0001, local_lTaylor=true)
            nphi = Q * ngrid
            phi_probe = build_sitewise_positive_state(setup.w_num, ngrid; seed=20260428)

            @testset "Stage 2: lifted periodic streaming operator satisfies the product rule" begin
                S_sparse = build_streaming_carleman_operator_sparse(Q, truncation_order, poly_order, ngrid; S_lbm=setup.S_lbm)
                block22, _ = carleman_C_block_dim(Q, 2, 2, 0)
                block33, _ = carleman_C_block_dim(Q, 3, 3, 0)

                Sphi = setup.S_lbm * phi_probe
                phi2 = state_kron(phi_probe, 2)
                phi3 = state_kron(phi_probe, 3)

                expected_order2 = kron(Sphi, phi_probe) + kron(phi_probe, Sphi)
                expected_order3 =
                    kron(Sphi, kron(phi_probe, phi_probe)) +
                    kron(phi_probe, kron(Sphi, phi_probe)) +
                    kron(phi_probe, kron(phi_probe, Sphi))

                actual_order2 = Array(S_sparse[block22, block22] * phi2)
                actual_order3 = Array(S_sparse[block33, block33] * phi3)

                @test actual_order2 ≈ expected_order2 atol=1e-11 rtol=1e-11
                @test actual_order3 ≈ expected_order3 atol=1e-11 rtol=1e-11

                println("✅ D2Q9 Stage 2 passed: lifted periodic streaming blocks reproduce the product-rule action on φ^[2] and φ^[3]")
            end

            @testset "Stage 3: streaming-only first level" begin
                V0 = Float64.(carleman_V(phi_probe, truncation_order))
                S_sparse = build_streaming_carleman_operator_sparse(Q, truncation_order, poly_order, ngrid; S_lbm=setup.S_lbm)

                direct_streaming_rhs = -setup.S_lbm * phi_probe
                lifted_streaming_rhs = Array(-S_sparse * V0)

                @test lifted_streaming_rhs[1:nphi] ≈ direct_streaming_rhs atol=1e-11 rtol=1e-11

                println("✅ D2Q9 Stage 3 passed: first Carleman level reproduces D2Q9 streaming-only RHS")
            end

            @testset "Stage 4: coupled collision + streaming first-level RHS" begin
                V0 = Float64.(carleman_V(phi_probe, truncation_order))
                C_sparse, bt_sparse, _ = carleman_C_sparse(Q, truncation_order, poly_order, nothing, nothing, tau_value, force_factor, w_value, e_value)
                S_sparse = build_streaming_carleman_operator_sparse(Q, truncation_order, poly_order, ngrid; S_lbm=setup.S_lbm)

                direct_rhs = direct_lbe_rhs_ngrid(phi_probe, setup.S_lbm, setup.F1_ngrid, setup.F2_ngrid, setup.F3_ngrid)
                lifted_rhs = Array((C_sparse - S_sparse) * V0 + bt_sparse)

                @test lifted_rhs[1:nphi] ≈ direct_rhs atol=1e-10 rtol=1e-10

                println("✅ D2Q9 Stage 4 passed: assembled coupled D2Q9 Carleman operator reproduces the direct first-level RHS")
            end

            @testset "Stage 5: one-step explicit Euler update" begin
                V0 = Float64.(carleman_V(phi_probe, truncation_order))
                C_sparse, bt_sparse, _ = carleman_C_sparse(Q, truncation_order, poly_order, nothing, nothing, tau_value, force_factor, w_value, e_value)
                S_sparse = build_streaming_carleman_operator_sparse(Q, truncation_order, poly_order, ngrid; S_lbm=setup.S_lbm)

                direct_rhs = direct_lbe_rhs_ngrid(phi_probe, setup.S_lbm, setup.F1_ngrid, setup.F2_ngrid, setup.F3_ngrid)
                lifted_rhs = Array((C_sparse - S_sparse) * V0 + bt_sparse)

                direct_next = phi_probe + dt * direct_rhs
                lifted_next = phi_probe + dt * lifted_rhs[1:nphi]

                @test lifted_next ≈ direct_next atol=1e-10 rtol=1e-10

                println("✅ D2Q9 Stage 5 passed: one explicit Euler step agrees between direct D2Q9 spatial LBE and CLBM first level")
            end
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

    if old_settings.F1_ngrid !== nothing
        global F1_ngrid = old_settings.F1_ngrid
    end
    if old_settings.F2_ngrid !== nothing
        global F2_ngrid = old_settings.F2_ngrid
    end
    if old_settings.F3_ngrid !== nothing
        global F3_ngrid = old_settings.F3_ngrid
    end
end

println("\n🎉 Streaming and collision+streaming Carleman validation completed successfully!")