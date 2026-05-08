using SymPy
using LinearAlgebra

function LBM_const_subs(omega, tau_value) 
    #=
    w1, w2, and w3 are weight factors for left moving particle f1, resting particle f0 (middle), and right moving particle, respectively  
    =#
    @syms tau

    w, e, w_value, e_value, a, b, c, d, a_value, b_value, c_value, d_value = lbm_const_sym()
    flat_e = e isa AbstractVector && !isempty(e) && first(e) isa AbstractVector ? collect(Iterators.flatten(e)) : collect(e)
    flat_e_value = e_value isa AbstractVector && !isempty(e_value) && first(e_value) isa AbstractVector ? collect(Iterators.flatten(e_value)) : collect(e_value)
    subs_vars = Any[tau]
    append!(subs_vars, flat_e)
    append!(subs_vars, collect(w))
    append!(subs_vars, [a, b, c, d])

    subs_values = Any[tau_value]
    append!(subs_values, flat_e_value)
    append!(subs_values, collect(w_value))
    append!(subs_values, [a_value, b_value, c_value, d_value])
    for i = 1:length(subs_vars)
        omega = omega.subs(subs_vars[i], subs_values[i])
    end
    return omega
end
