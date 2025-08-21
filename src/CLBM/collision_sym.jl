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
    u = sum(e .* f) / rho
    return rho, u
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

    if lTaylor == true
        u = sum(e .* f) * (2-rho)
    else
        u = sum(e .* f) / rho
    end
    
    u = expand(u)
    u = collect(u,f)

    if lTaylor == true
        sum_e_f = sum(e .* f)
        eiu = e .* sum_e_f
        eiu2 = eiu .^2
      #  feq = w .* (rho .+ (3/c) .* eiu + (9/(2*c^2)) .* (2-rho) .* (eiu2) .- (3/(2*c^2)) .* (2-rho) .* sum_e_f .^2) 
        feq = w .* (a * rho .+ b .* eiu + c .* (2-rho) .* (eiu2) .+ d .* (2-rho) .* sum_e_f .^2) 
    else
        eiu = e.*u
        eiu2 = eiu.^2
        u2 = sum(u.*u)
       # feq = w.*rho .+ rho.*(w.*((3/c).*eiu+(9/(2*c^2)).*(eiu2).-(3/(2*c^2)).*u2)) 
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

