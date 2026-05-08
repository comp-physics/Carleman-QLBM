using LinearAlgebra
using SparseArrays

# Semi-discrete n-point LBE (Li et al., Eq. eq:eq:lbm-carleman-n)
#
#   ∂_t φ(x)  =  −S φ(x)  +  F^(1) φ(x)
#                         +  F^(2) φ^[2](x)
#                         +  F^(3) φ^[3](x)
#
# Streaming `S` is a centered-difference operator (or whatever finite-difference
# stencil the caller chooses). The collision is the polynomial Taylor-expanded
# form obtained with `lTaylor = true` (so that 1/ρ ≈ 2 − ρ keeps the operator
# polynomial in φ).
#
# This module is intentionally independent of the Carleman lift: it is the
# direct nonlinear ODE that CLBE linearizes. Its role in the repository:
#   * primary self-consistency target for CLBE (the apples-to-apples reference),
#   * bridge between the discrete LBM baseline (src/LBM/, exact lattice shift,
#     integer iteration) and the Carleman-linearized system (src/CLBE/).
#
# Inputs are passed by caller — no globals, no fallback streaming. This keeps
# the module standalone and composable with any streaming stencil supplied
# elsewhere in the codebase.

function direct_lbe_rhs_ngrid(phi, S_lbm, F1_ngrid, F2_ngrid, F3_ngrid)
    rhs = -S_lbm * phi + F1_ngrid * phi
    if F2_ngrid !== nothing
        rhs += F2_ngrid * kron(phi, phi)
    end
    if F3_ngrid !== nothing
        rhs += F3_ngrid * kron(phi, kron(phi, phi))
    end
    return rhs
end

function normalize_direct_lbe_integrator(integrator)
    if integrator in (:euler, "euler")
        return :euler
    elseif integrator in (:exponential_euler, :etd, :etd1, "exponential_euler", "etd", "etd1")
        return :exponential_euler
    else
        error("Unsupported direct LBE integrator $(repr(integrator)). Supported options are :euler and :exponential_euler (aliases :etd, :etd1).")
    end
end

function direct_lbe_linear_mul(phi, S_lbm, F1_ngrid)
    return -S_lbm * phi + F1_ngrid * phi
end

function direct_lbe_nonlinear_rhs_ngrid(phi, F2_ngrid, F3_ngrid)
    rhs = zeros(Float64, length(phi))
    if F2_ngrid !== nothing
        rhs .+= F2_ngrid * kron(phi, phi)
    end
    if F3_ngrid !== nothing
        rhs .+= F3_ngrid * kron(phi, kron(phi, phi))
    end
    return rhs
end

function direct_lbe_affine_augmented_mul(S_lbm, F1_ngrid, nonlinear_term, x)
    state_dim = length(nonlinear_term)
    y = zeros(Float64, state_dim + 1)
    y[1:state_dim] .= direct_lbe_linear_mul(view(x, 1:state_dim), S_lbm, F1_ngrid) .+ nonlinear_term .* x[end]
    y[end] = 0.0
    return y
end

function krylov_expv_direct_lbe_affine(S_lbm, F1_ngrid, nonlinear_term, v, dt; m=30, tol=1e-10)
    n_aug = length(v)
    beta = norm(v)
    if iszero(beta)
        return copy(v)
    end

    m_eff = min(m, n_aug)
    V = zeros(Float64, n_aug, m_eff + 1)
    H = zeros(Float64, m_eff + 1, m_eff)
    V[:, 1] = v ./ beta

    actual_m = m_eff
    for j = 1:m_eff
        w = direct_lbe_affine_augmented_mul(S_lbm, F1_ngrid, nonlinear_term, view(V, :, j))
        for i = 1:j
            hij = dot(view(V, :, i), w)
            H[i, j] = hij
            w .-= hij .* view(V, :, i)
        end

        hnext = norm(w)
        H[j + 1, j] = hnext
        if hnext <= tol
            actual_m = j
            break
        end

        if j < m_eff
            V[:, j + 1] = w ./ hnext
        end
    end

    H_small = H[1:actual_m, 1:actual_m]
    e1 = zeros(Float64, actual_m)
    e1[1] = beta
    y_small = exp(dt * H_small) * e1
    return V[:, 1:actual_m] * y_small
end

function direct_lbe_exponential_euler_step(phi, dt, S_lbm, F1_ngrid, F2_ngrid, F3_ngrid)
    nonlinear_term = direct_lbe_nonlinear_rhs_ngrid(phi, F2_ngrid, F3_ngrid)
    augmented_state = vcat(Float64.(phi), 1.0)
    propagated = krylov_expv_direct_lbe_affine(S_lbm, F1_ngrid, nonlinear_term, augmented_state, dt)
    return propagated[1:length(phi)]
end

function timeMarching_direct_LBE_ngrid(phi_ini, dt, n_time, F1_ngrid, F2_ngrid, F3_ngrid; S_lbm, integrator=:euler)
    integrator_key = normalize_direct_lbe_integrator(integrator)
    phiT = zeros(length(phi_ini), n_time)
    phiT[:, 1] = Float64.(phi_ini)

    for nt = 2:n_time
        if integrator_key == :euler
            rhs = direct_lbe_rhs_ngrid(phiT[:, nt - 1], S_lbm, F1_ngrid, F2_ngrid, F3_ngrid)
            phiT[:, nt] = phiT[:, nt - 1] + dt * rhs
        else
            phiT[:, nt] = direct_lbe_exponential_euler_step(phiT[:, nt - 1], dt, S_lbm, F1_ngrid, F2_ngrid, F3_ngrid)
        end
    end

    return phiT
end
