using SymPy
using LinearAlgebra

function LBM_const_subs(omega, tau_value) 
    #=
    w1, w2, and w3 are weight factors for left moving particle f1, resting particle f0 (middle), and right moving particle, respectively  
    =#
    @syms tau

    w, e, w_value, e_value, a, b, c, d, a_value, b_value, c_value, d_value = lbm_const_sym()
    subs_vars = collect(Iterators.flatten([tau, e, w, a, b, c, d]))
    subs_values = collect(Iterators.flatten([tau_value, e_value, w_value, a_value, b_value, c_value, d_value])) 
    for i = 1:length(subs_vars)
        omega = omega.subs(subs_vars[i], subs_values[i])
    end
    return omega
end
