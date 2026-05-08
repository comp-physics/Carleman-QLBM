if !haskey(ENV, "QCFD_HOME")
    ENV["QCFD_HOME"] = normpath(joinpath(@__DIR__, "..", ".."))
end
if !haskey(ENV, "QCFD_SRC")
    ENV["QCFD_SRC"] = joinpath(ENV["QCFD_HOME"], "src") * "/"
end

QCFD_SRC = ENV["QCFD_SRC"]

using LinearAlgebra
using SparseArrays

include(QCFD_SRC * "LBM/lbm_const_sym.jl")
include(QCFD_SRC * "LBM/cal_feq.jl")

function index2(i, j, Q)
    return (i - 1) * Q + j
end

function index3(i, j, k, Q)
    return (i - 1) * Q^2 + (j - 1) * Q + k
end

function solve_vector_coefficients(system_matrix, rhs_vectors)
    rhs_matrix = reduce(vcat, permutedims.(rhs_vectors))
    coeff_matrix = system_matrix \ rhs_matrix
    return [vec(coeff_matrix[row, :]) for row = 1:size(coeff_matrix, 1)]
end

function lbm_const_numerical(; Q_local, D_local)
    if Q_local == 3 && D_local == 1
        w_value = [1.0 / 6, 2.0 / 3, 1.0 / 6]
        e_value = [-1.0, 0.0, 1.0]
    elseif Q_local == 9 && D_local == 2
        w_value = [4.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 36, 1.0 / 36, 1.0 / 36, 1.0 / 36]
        e_value = [
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
    else
        error("Unsupported lattice for numerical direct-LBE coefficient generation: Q = $Q_local, D = $D_local")
    end

    a_value, b_value, c_value, d_value = 1.0, 3.0, 9.0 / 2.0, -3.0 / 2.0
    return w_value, e_value, a_value, b_value, c_value, d_value
end

function numerical_direct_lbe_collision_rhs(state, tau_value; w_value_input, e_value_input, rho_value_input=1.0, lTaylor_input=true, D_input)
    Q_local = length(state)
    D_local = D_input
    w_local, e_local_default, a_value, b_value, c_value, d_value = lbm_const_numerical(Q_local=Q_local, D_local=D_local)
    w_local = Float64.(w_value_input === nothing ? w_local : w_value_input)
    e_local = e_value_input === nothing ? e_local_default : e_value_input

    if !lTaylor_input && rho_value_input != 1.0
        error("Numerical direct-LBE coefficient generation requires lTaylor=true or rho_value=1.0 so the collision operator remains polynomial.")
    end

    rho_local = rho_value_input == 1.0 ? rho_value_input : sum(state)
    momentum = D_local == 1 ? sum(e_local .* state) : [sum(ei[dim] * fi for (ei, fi) in zip(e_local, state)) for dim = 1:D_local]

    if lTaylor_input
        sum_e_f = momentum
        eiu = D_local == 1 ? e_local .* sum_e_f : [sum(ei[d] * sum_e_f[d] for d = 1:D_local) for ei in e_local]
        eiu2 = eiu .^ 2
        momentum_sq = D_local == 1 ? sum_e_f^2 : sum(component^2 for component in sum_e_f)
        feq = w_local .* (a_value * rho_local .+ b_value .* eiu + c_value .* (2.0 - rho_local) .* eiu2 .+ d_value .* (2.0 - rho_local) .* momentum_sq)
    else
        velocity = D_local == 1 ? momentum / rho_local : [momentum[dim] / rho_local for dim = 1:D_local]
        eiu = D_local == 1 ? e_local .* velocity : [sum(ei[d] * velocity[d] for d = 1:D_local) for ei in e_local]
        eiu2 = eiu .^ 2
        velocity_sq = D_local == 1 ? velocity^2 : sum(component^2 for component in velocity)
        feq = w_local .* rho_local .+ rho_local .* (w_local .* (b_value .* eiu .+ c_value .* eiu2 .+ d_value .* velocity_sq))
    end

    return Float64.((feq .- state) ./ tau_value)
end

function numerical_direct_lbe_coefficients(poly_order, Q, tau_value; w_value_input, e_value_input, rho_value_input=1.0, lTaylor_input=true, D_input)
    F1 = zeros(Q, Q)
    F2 = poly_order >= 2 ? zeros(Q, Q^2) : nothing
    F3 = poly_order >= 3 ? zeros(Q, Q^3) : nothing

    evaluate_state(entries) = begin
        state = zeros(Q)
        for (ind, value) in entries
            state[ind] = value
        end
        numerical_direct_lbe_collision_rhs(
            state,
            tau_value;
            w_value_input=w_value_input,
            e_value_input=e_value_input,
            rho_value_input=rho_value_input,
            lTaylor_input=lTaylor_input,
            D_input=D_input,
        )
    end

    one_var_system = [1.0 1.0 1.0; 2.0 4.0 8.0; 3.0 9.0 27.0]
    linear_coeffs = zeros(Q, Q)
    quadratic_diagonal = poly_order >= 2 ? zeros(Q, Q) : nothing
    cubic_diagonal = poly_order >= 3 ? zeros(Q, Q) : nothing
    quadratic_cross = poly_order >= 2 ? Dict{Tuple{Int,Int}, Vector{Float64}}() : nothing
    cubic_iij = poly_order >= 3 ? Dict{Tuple{Int,Int}, Vector{Float64}}() : nothing
    cubic_ijj = poly_order >= 3 ? Dict{Tuple{Int,Int}, Vector{Float64}}() : nothing

    for i = 1:Q
        rhs_vectors = [
            evaluate_state([(i, 1.0)]),
            evaluate_state([(i, 2.0)]),
            evaluate_state([(i, 3.0)]),
        ]
        coeffs = solve_vector_coefficients(one_var_system, rhs_vectors)
        linear_coeffs[:, i] .= coeffs[1]
        F1[:, i] .= coeffs[1]

        if poly_order >= 2
            quadratic_diagonal[:, i] .= coeffs[2]
            F2[:, index2(i, i, Q)] .= coeffs[2]
        end

        if poly_order >= 3
            cubic_diagonal[:, i] .= coeffs[3]
            F3[:, index3(i, i, i, Q)] .= coeffs[3]
        end
    end

    if poly_order >= 2
        pair_system = [1.0 1.0 1.0; 2.0 4.0 2.0; 2.0 2.0 4.0]
        for i = 1:Q-1
            for j = i+1:Q
                residual(x, y) = begin
                    value = evaluate_state([(i, x), (j, y)])
                    value .-= linear_coeffs[:, i] .* x
                    value .-= linear_coeffs[:, j] .* y
                    value .-= quadratic_diagonal[:, i] .* x^2
                    value .-= quadratic_diagonal[:, j] .* y^2
                    if poly_order >= 3
                        value .-= cubic_diagonal[:, i] .* x^3
                        value .-= cubic_diagonal[:, j] .* y^3
                    end
                    value
                end

                coeffs = solve_vector_coefficients(pair_system, [
                    residual(1.0, 1.0),
                    residual(2.0, 1.0),
                    residual(1.0, 2.0),
                ])

                quadratic_cross[(i, j)] = coeffs[1]
                F2[:, index2(i, j, Q)] .= coeffs[1]

                if poly_order >= 3
                    cubic_iij[(i, j)] = coeffs[2]
                    cubic_ijj[(i, j)] = coeffs[3]
                    F3[:, index3(i, i, j, Q)] .= coeffs[2]
                    F3[:, index3(i, j, j, Q)] .= coeffs[3]
                end
            end
        end
    end

    if poly_order >= 3
        for i = 1:Q-2
            for j = i+1:Q-1
                for k = j+1:Q
                    residual = evaluate_state([(i, 1.0), (j, 1.0), (k, 1.0)])
                    residual .-= linear_coeffs[:, i] .+ linear_coeffs[:, j] .+ linear_coeffs[:, k]
                    residual .-= quadratic_diagonal[:, i] .+ quadratic_diagonal[:, j] .+ quadratic_diagonal[:, k]
                    residual .-= cubic_diagonal[:, i] .+ cubic_diagonal[:, j] .+ cubic_diagonal[:, k]

                    for pair in ((i, j), (i, k), (j, k))
                        residual .-= quadratic_cross[pair]
                        residual .-= cubic_iij[pair]
                        residual .-= cubic_ijj[pair]
                    end

                    F3[:, index3(i, j, k, Q)] .= residual
                end
            end
        end
    end

    return F1, F2, F3
end

function sparse_lift_local_collision(F_single::AbstractMatrix, Q::Int, ngrid::Int, order::Int)
    nQ = ngrid * Q
    row_dim = nQ
    col_dim = nQ^order
    I = Int[]
    J = Int[]
    V = Float64[]

    if order == 1
        for α = 1:ngrid, r = 1:Q, c = 1:Q
            val = F_single[r, c]
            iszero(val) && continue
            push!(I, (α - 1) * Q + r)
            push!(J, (α - 1) * Q + c)
            push!(V, val)
        end
    elseif order == 2
        for α = 1:ngrid, r = 1:Q, c1 = 1:Q, c2 = 1:Q
            val = F_single[r, index2(c1, c2, Q)]
            iszero(val) && continue
            g1 = (α - 1) * Q + c1
            g2 = (α - 1) * Q + c2
            push!(I, (α - 1) * Q + r)
            push!(J, index2(g1, g2, nQ))
            push!(V, val)
        end
    elseif order == 3
        for α = 1:ngrid, r = 1:Q, c1 = 1:Q, c2 = 1:Q, c3 = 1:Q
            val = F_single[r, index3(c1, c2, c3, Q)]
            iszero(val) && continue
            g1 = (α - 1) * Q + c1
            g2 = (α - 1) * Q + c2
            g3 = (α - 1) * Q + c3
            push!(I, (α - 1) * Q + r)
            push!(J, index3(g1, g2, g3, nQ))
            push!(V, val)
        end
    else
        error("Unsupported order = $order")
    end

    return sparse(I, J, V, row_dim, col_dim)
end

function build_sparse_ngrid_coefficients(poly_order, Q, tau_value, ngrid; w_value_input, e_value_input, rho_value_input=1.0, lTaylor_input=true, D_input)
    F1_single, F2_single, F3_single = numerical_direct_lbe_coefficients(
        poly_order,
        Q,
        tau_value;
        w_value_input=w_value_input,
        e_value_input=e_value_input,
        rho_value_input=rho_value_input,
        lTaylor_input=lTaylor_input,
        D_input=D_input,
    )

    F1_ngrid_local = sparse_lift_local_collision(F1_single, Q, ngrid, 1)
    F2_ngrid_local = poly_order >= 2 ? sparse_lift_local_collision(F2_single, Q, ngrid, 2) : nothing
    F3_ngrid_local = poly_order >= 3 ? sparse_lift_local_collision(F3_single, Q, ngrid, 3) : nothing
    return F1_ngrid_local, F2_ngrid_local, F3_ngrid_local
end

function streaming_operator_D2Q9_interleaved_periodic(nx, ny, hx, hy)
    e = [
        [0, 0],
        [1, 0],
        [-1, 0],
        [0, 1],
        [0, -1],
        [1, 1],
        [1, -1],
        [-1, 1],
        [-1, -1],
    ]

    n_velocities = 9
    n_total = n_velocities * nx * ny
    I_idx = Int[]
    J_idx = Int[]
    V_vals = Float64[]
    global_index(vel_idx, i, j) = ((j - 1) * nx + (i - 1)) * n_velocities + vel_idx

    for j in 1:ny
        for i in 1:nx
            for vel in 1:n_velocities
                row_idx = global_index(vel, i, j)
                ex, ey = e[vel]

                if ex == 0 && ey == 0
                    push!(I_idx, row_idx); push!(J_idx, row_idx); push!(V_vals, 0.0)
                    continue
                end

                if ex != 0
                    left_i = i == 1 ? nx : i - 1
                    right_i = i == nx ? 1 : i + 1
                    push!(I_idx, row_idx); push!(J_idx, global_index(vel, left_i, j)); push!(V_vals, -ex / (2 * hx))
                    push!(I_idx, row_idx); push!(J_idx, global_index(vel, right_i, j)); push!(V_vals, ex / (2 * hx))
                end

                if ey != 0
                    bottom_j = j == 1 ? ny : j - 1
                    top_j = j == ny ? 1 : j + 1
                    push!(I_idx, row_idx); push!(J_idx, global_index(vel, i, bottom_j)); push!(V_vals, -ey / (2 * hy))
                    push!(I_idx, row_idx); push!(J_idx, global_index(vel, i, top_j)); push!(V_vals, ey / (2 * hy))
                end
            end
        end
    end

    return sparse(I_idx, J_idx, V_vals, n_total, n_total), e
end

function opposite_velocity_index_D2Q9(vel)
    opposite_map = Dict(1 => 1, 2 => 3, 3 => 2, 4 => 5, 5 => 4, 6 => 9, 7 => 8, 8 => 7, 9 => 6)
    return opposite_map[vel]
end

function streaming_operator_D2Q9_interleaved_boundary_aware(nx, ny, hx, hy)
    e = [
        [0, 0],
        [1, 0],
        [-1, 0],
        [0, 1],
        [0, -1],
        [1, 1],
        [1, -1],
        [-1, 1],
        [-1, -1],
    ]

    n_velocities = 9
    n_total = n_velocities * nx * ny
    I_idx = Int[]
    J_idx = Int[]
    V_vals = Float64[]
    global_index(vel_idx, i, j) = ((j - 1) * nx + (i - 1)) * n_velocities + vel_idx

    for j in 1:ny
        for i in 1:nx
            for vel in 1:n_velocities
                row_idx = global_index(vel, i, j)
                ex, ey = e[vel]

                if ex == 0 && ey == 0
                    push!(I_idx, row_idx); push!(J_idx, row_idx); push!(V_vals, 0.0)
                    continue
                end

                if ex != 0
                    left_i = i == 1 ? nx : i - 1
                    right_i = i == nx ? 1 : i + 1
                    push!(I_idx, row_idx); push!(J_idx, global_index(vel, left_i, j)); push!(V_vals, -ex / (2 * hx))
                    push!(I_idx, row_idx); push!(J_idx, global_index(vel, right_i, j)); push!(V_vals, ex / (2 * hx))
                end

                if ey != 0
                    reflected_vel = opposite_velocity_index_D2Q9(vel)
                    if j > 1
                        push!(I_idx, row_idx); push!(J_idx, global_index(vel, i, j - 1)); push!(V_vals, -ey / (2 * hy))
                    else
                        push!(I_idx, row_idx); push!(J_idx, global_index(reflected_vel, i, j)); push!(V_vals, -ey / (2 * hy))
                    end

                    if j < ny
                        push!(I_idx, row_idx); push!(J_idx, global_index(vel, i, j + 1)); push!(V_vals, ey / (2 * hy))
                    else
                        push!(I_idx, row_idx); push!(J_idx, global_index(reflected_vel, i, j)); push!(V_vals, ey / (2 * hy))
                    end
                end
            end
        end
    end

    return sparse(I_idx, J_idx, V_vals, n_total, n_total), e
end

function select_d2q9_direct_lbe_streaming_operator(nx, ny, hx, hy; boundary_setup=false)
    if boundary_setup
        return streaming_operator_D2Q9_interleaved_boundary_aware(nx, ny, hx, hy)
    end
    return streaming_operator_D2Q9_interleaved_periodic(nx, ny, hx, hy)
end

function tg_velocity_field(i, j, nx, ny, amplitude)
    x = 2 * π * (i - 1) / nx
    y = 2 * π * (j - 1) / ny
    ux = amplitude * sin(x) * cos(y)
    uy = -amplitude * cos(x) * sin(y)
    return ux, uy
end

function d2q9_equilibrium(rho, ux, uy, w_value, e_value, a_value, b_value, c_value, d_value)
    vcx = [velocity[1] for velocity in e_value]
    vcy = [velocity[2] for velocity in e_value]
    return cal_feq(rho, ux, uy, w_value, vcx, vcy, a_value, b_value, c_value, d_value)
end

function tg2d_initial_condition(nx, ny, amplitude, rho_value, w_value, e_value, a_value, b_value, c_value, d_value)
    phi_ini = zeros(9 * nx * ny)
    for j in 1:ny
        for i in 1:nx
            ux, uy = tg_velocity_field(i, j, nx, ny, amplitude)
            feq = d2q9_equilibrium(rho_value, ux, uy, w_value, e_value, a_value, b_value, c_value, d_value)
            start_idx = ((j - 1) * nx + (i - 1)) * 9 + 1
            phi_ini[start_idx:start_idx + 8] .= feq
        end
    end
    return phi_ini
end

function macroscopic_fields_from_state(phi, nx, ny, e_value)
    rho = zeros(nx, ny)
    ux = zeros(nx, ny)
    uy = zeros(nx, ny)

    for j in 1:ny
        for i in 1:nx
            start_idx = ((j - 1) * nx + (i - 1)) * 9 + 1
            fi = phi[start_idx:start_idx + 8]
            rho_local = sum(fi)
            rho[i, j] = rho_local
            ux[i, j] = sum(fi[k] * e_value[k][1] for k = 1:9) / rho_local
            uy[i, j] = sum(fi[k] * e_value[k][2] for k = 1:9) / rho_local
        end
    end

    return rho, ux, uy
end

function build_direct_lbe_d2q9_runtime(; nx, ny, rho_value=1.0001, tau_value=1.0, boundary_setup=false, hx=1.0, hy=1.0, poly_order_input=3, lTaylor_input=true)
    Q_local = 9
    D_local = 2
    nspatial = nx * ny

    _, _, numeric_weights, numeric_velocities, _, _, _, _, a_val, b_val, c_val, d_val = lbm_const_sym(Q_local=Q_local, D_local=D_local)
    F1_ngrid, F2_ngrid, F3_ngrid = build_sparse_ngrid_coefficients(
        poly_order_input,
        Q_local,
        tau_value,
        nspatial;
        w_value_input=numeric_weights,
        e_value_input=numeric_velocities,
        rho_value_input=rho_value,
        lTaylor_input=lTaylor_input,
        D_input=D_local,
    )

    S_lbm, e_stream = select_d2q9_direct_lbe_streaming_operator(nx, ny, hx, hy; boundary_setup=boundary_setup)

    return (
        Q=Q_local,
        D=D_local,
        ngrid=nspatial,
        numeric_weights=numeric_weights,
        numeric_velocities=numeric_velocities,
        a_val=a_val,
        b_val=b_val,
        c_val=c_val,
        d_val=d_val,
        F1_ngrid=F1_ngrid,
        F2_ngrid=F2_ngrid,
        F3_ngrid=F3_ngrid,
        S_lbm=S_lbm,
        e_value=numeric_velocities,
        streaming_velocities=e_stream,
    )
end

function build_direct_lbe_tg_reference(; nx, ny, amplitude, rho_value, local_n_time, boundary_setup, runtime, dt_value, direct_lbe_integrator=:euler)
    phi_ini = tg2d_initial_condition(
        nx,
        ny,
        amplitude,
        rho_value,
        runtime.numeric_weights,
        runtime.e_value,
        runtime.a_val,
        runtime.b_val,
        runtime.c_val,
        runtime.d_val,
    )

    reference_phi_history = timeMarching_direct_LBE_ngrid(
        phi_ini,
        dt_value,
        local_n_time,
        runtime.F1_ngrid,
        runtime.F2_ngrid,
        runtime.F3_ngrid;
        S_lbm=runtime.S_lbm,
        integrator=direct_lbe_integrator,
    )

    return (
        reference_initial_state=phi_ini,
        reference_phi_history=reference_phi_history,
    )
end