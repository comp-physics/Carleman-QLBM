l_sympy = true

QCFD_HOME = ENV["QCFD_HOME"]  
include(QCFD_HOME * "/julia_lib/matrix_kit.jl")

if l_sympy
    using SymPy
else
    using Symbolics
end

function lbm_u(e, f)
    rho = sum(f)
    if e isa AbstractVector && !isempty(e) && first(e) isa AbstractVector
        u = [sum(ei[d] * fi for (ei, fi) in zip(e, f)) / rho for d = 1:length(first(e))]
    else
        u = sum(e .* f) / rho
    end
    return rho, u
end

function lbm_dot_velocity(e, vector_state)
    return [sum(ei[d] * vector_state[d] for d = 1:length(vector_state)) for ei in e]
end

function lbm_expand_state(state, f)
    if state isa AbstractVector
        return [collect(expand(component), f) for component in state]
    end
    return collect(expand(state), f)
end

function collision(Q, D, w, e, rho, lTaylor, lorder2)
    #Declare your discrete density variables
    if l_sympy
        f = [symbols("f$i") for i in 1:Q]
        tau = symbols("tau");
    else
        f =  Symbolics.variables(:f, 1:Q)
       # @variables tau
        tau = Symbolics.variable("tau")
    end

    #---LBM constant---
_, _, _, _, a, b, c, d, a_value, b_value, c_value, d_value = lbm_const_sym()
    #

    #Assume incompressible flow with unity constant density
    if rho!=1
       rho = sum(f); #2022-09-21/XYLI
    end

    momentum = if D == 1
        sum(e .* f)
    else
        [sum(ei[dim] * fi for (ei, fi) in zip(e, f)) for dim = 1:D]
    end

    if lTaylor == true
        u = D == 1 ? momentum * (2 - rho) : [momentum[dim] * (2 - rho) for dim = 1:D]
    else
        u = D == 1 ? momentum / rho : [momentum[dim] / rho for dim = 1:D]
    end

    u = lbm_expand_state(u, f)

    if lTaylor == true
        sum_e_f = momentum
        eiu = D == 1 ? e .* sum_e_f : lbm_dot_velocity(e, sum_e_f)
        eiu2 = eiu .^2
        momentum_sq = D == 1 ? sum_e_f .^ 2 : sum(component^2 for component in sum_e_f)
        feq = w .* (a * rho .+ b .* eiu + c .* (2-rho) .* (eiu2) .+ d .* (2-rho) .* momentum_sq)
    else
        eiu = D == 1 ? e .* u : lbm_dot_velocity(e, u)
        eiu2 = eiu.^2
        u2 = D == 1 ? u^2 : sum(component^2 for component in u)
        feq = w .* rho .+ rho .* (w .* (b .* eiu .+ c .* (eiu2) .+ d .* u2))
    end
    feq = [expand(i) for i in feq]

    #Calculate the collision term driving the differential equation
#    omega = -(dt ./ tau) .* (f .- feq)
    omega = -(1 ./ tau) .* (f .- feq)
#    if l_sympy
#        omega_sub = omega.subs(tau, 1)
#    else
#        omega_sub = substitute(omega, Dict(tau=>1.))
#    end
    #
    return f, omega, u, rho
end

