QCFD_SRC = ENV["QCFD_SRC"]  
QCFD_HOME = ENV["QCFD_HOME"]  

using SparseArrays

include(QCFD_SRC * "CLBE/carleman_transferA.jl")

function infer_lbm_dimension(e_value)
    if e_value isa AbstractVector && !isempty(e_value) && first(e_value) isa AbstractVector
        return length(first(e_value))
    end
    return 1
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
        error("Unsupported lattice for numerical coefficient generation: Q = $Q_local, D = $D_local")
    end

    a_value, b_value, c_value, d_value = 1.0, 3.0, 9.0 / 2.0, -3.0 / 2.0
    return w_value, e_value, a_value, b_value, c_value, d_value
end

function numerical_collision_rhs(state, tau_value; w_value_input=nothing, e_value_input=nothing, rho_value_input=nothing, lTaylor_input=nothing, D_input=nothing)
    Q_local = length(state)
    e_local = e_value_input === nothing ? (@isdefined(e_value) ? e_value : nothing) : e_value_input
    D_local = D_input === nothing ? (e_local === nothing ? (@isdefined(D) ? D : 1) : infer_lbm_dimension(e_local)) : D_input

    w_local, e_default, a_value, b_value, c_value, d_value = lbm_const_numerical(Q_local=Q_local, D_local=D_local)
    w_local = w_value_input === nothing ? w_local : Float64.(w_value_input)
    e_local = e_local === nothing ? e_default : e_local

    rho_reference = rho_value_input === nothing ? (@isdefined(rho0) ? rho0 : 1.0) : rho_value_input
    lTaylor_local = lTaylor_input === nothing ? (@isdefined(lTaylor) ? lTaylor : false) : lTaylor_input

    if !lTaylor_local && rho_reference != 1.0
        error("Numerical coefficient generation currently requires lTaylor=true or rho_value=1.0 so the collision operator remains polynomial.")
    end

    rho_local = rho_reference == 1.0 ? rho_reference : sum(state)

    momentum = if D_local == 1
        sum(e_local .* state)
    else
        [sum(ei[dim] * fi for (ei, fi) in zip(e_local, state)) for dim = 1:D_local]
    end

    if lTaylor_local
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

function numerical_carleman_coefficients(poly_order, Q, tau_value; w_value_input=nothing, e_value_input=nothing, rho_value_input=nothing, lTaylor_input=nothing, D_input=nothing)
    F1 = zeros(Q, Q)
    F2 = poly_order >= 2 ? zeros(Q, Q^2) : nothing
    F3 = poly_order >= 3 ? zeros(Q, Q^3) : nothing

    evaluate_state(entries) = begin
        state = zeros(Q)
        for (ind, value) in entries
            state[ind] = value
        end
        numerical_collision_rhs(
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


function cal_delta_alpha(which_alpha, ngrid)
    delta_alpha = zeros(ngrid)
    delta_alpha[which_alpha] = 1.0
    return delta_alpha
end

function local_selector_matrix(which_alpha, Q, ngrid)
    selector = zeros(Q, Q * ngrid)
    ind_s = (which_alpha - 1) * Q + 1
    ind_e = which_alpha * Q
    selector[:, ind_s:ind_e] = Matrix(1.0I, Q, Q)
    return selector
end

function coeff_LBM_Fi_xalpha(Fj, which_alpha, ngrid, order)
    # F^(k)(x_alpha) acts on phi^[k] with dimension (nQ)^k, while remaining local in x.
    selector_alpha = local_selector_matrix(which_alpha, size(Fj, 1), ngrid)
    selector_alpha_k = Kron_kth(selector_alpha, order)
    F_xalpha = Fj * selector_alpha_k
    return F_xalpha
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

function build_sparse_ngrid_coefficients(poly_order, Q, tau_value, ngrid; w_value_input=nothing, e_value_input=nothing, rho_value_input=nothing, lTaylor_input=nothing, D_input=nothing)
    F1_single, F2_single, F3_single = numerical_carleman_coefficients(
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

function coeff_LBM_Fi_ngrid(Q, j, f, omega, tau_value, ngrid; method=:symbolic, w_value_input=nothing, e_value_input=nothing, rho_value_input=nothing, lTaylor_input=nothing, D_input=nothing)
    row_dim = ngrid * Q
    col_dim = (ngrid * Q)^j
    Fj_ngrid = zeros(row_dim, col_dim)
    # F^(k)(x_alpha) is a Q × (nQ)^k matrix selecting only local monomials at x_alpha.
    # F^(k)(x) is the n-point collision operator with size nQ × (nQ)^k.
    if method == :symbolic
        Fj = F_carlemanOrder_collision(Q, j, f, omega, tau_value)
    elseif method == :numerical
        F1_single, F2_single, F3_single = numerical_carleman_coefficients(
            j,
            Q,
            tau_value;
            w_value_input=w_value_input,
            e_value_input=e_value_input,
            rho_value_input=rho_value_input,
            lTaylor_input=lTaylor_input,
            D_input=D_input,
        )
        Fj = j == 1 ? F1_single : (j == 2 ? F2_single : F3_single)
        return Fj, sparse_lift_local_collision(Fj, Q, ngrid, j)
    else
        error("Unknown coefficient-generation method: $method")
    end
    for i = 1 : ngrid
        F_xalpha = coeff_LBM_Fi_xalpha(Fj, i, ngrid, j)
        ind_s = (i - 1) * Q + 1
        ind_e = i * Q
        Fj_ngrid[ind_s:ind_e,:] = F_xalpha
    end
    return Fj, Fj_ngrid
end

function get_coeff_LBM_Fi_ngrid(poly_order, Q, f, omega, tau_value, ngrid; method=:symbolic, w_value_input=nothing, e_value_input=nothing, rho_value_input=nothing, lTaylor_input=nothing, D_input=nothing)
    if poly_order < 1 || poly_order > 3
        error("poly_order must be 1, 2, or 3")
    end

    if method == :numerical
        return build_sparse_ngrid_coefficients(
            poly_order,
            Q,
            tau_value,
            ngrid;
            w_value_input=w_value_input,
            e_value_input=e_value_input,
            rho_value_input=rho_value_input,
            lTaylor_input=lTaylor_input,
            D_input=D_input,
        )
    end

    _, F1_ngrid = coeff_LBM_Fi_ngrid(Q, 1, f, omega, tau_value, ngrid; method=method)
    F2_ngrid = nothing
    F3_ngrid = nothing

    if poly_order >= 2
        _, F2_ngrid = coeff_LBM_Fi_ngrid(Q, 2, f, omega, tau_value, ngrid; method=method)
    end

    if poly_order >= 3
        _, F3_ngrid = coeff_LBM_Fi_ngrid(Q, 3, f, omega, tau_value, ngrid; method=method)
    end

    return F1_ngrid, F2_ngrid, F3_ngrid
end

function transferA_ngrid(i, j, Q, ngrid)
    # dim(A_ij) = (ngrid*Q)^(i-1) * Q^j
    if j == 1
        Fj_ngrid = F1_ngrid
    elseif j == 2
        Fj_ngrid = F2_ngrid
    elseif j == 3
        Fj_ngrid = F3_ngrid
    else
        error("j of F^{j} must be 1, 2, 3, ..., poly_order")
    end
    #
    A_ij = sum_Kron_kth_identity(Fj_ngrid, i, Q * ngrid)
    #
    return A_ij
end

function transferA_S(i, Q, ngrid, S_Fj)
    #
    A_ij = sum_Kron_kth_identity(S_Fj, i, Q * ngrid)
    #
    return A_ij
end


function get_phi(f, ngrid)
    return vcat([collect(f) for _ = 1:ngrid]...)
end

function get_S_Fj(S, ngrid)
    # Streaming is already the full nonlocal n-point linear operator S on phi(x).
    return S
end
