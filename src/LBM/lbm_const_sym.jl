using SymPy

function lbm_const_sym()
    #=
    w1, w2, and w3 are weight factors for left moving particle f1, resting particle f0 (middle), and right moving particle, respectively  
    =#
    @syms e1, e2, e3, w1, w2, w3, a, b, c, d

    w = [w1, w2, w3]; #D1Q3
    w_value = [1. /6, 2. /3, 1. /6]; #D1Q3

    e = [e1, e2, e3]; #D1Q3
    e_value = [-1.0, 0.0, 1.0]; #D1Q3

    a_value, b_value, c_value, d_value = 1., 3., 9. / 2., - 3. / 2.

    return w, e, w_value, e_value, a, b, c, d, a_value, b_value, c_value, d_value
end
