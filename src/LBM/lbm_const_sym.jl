using SymPy

function lbm_const_sym(; Q_local=nothing, D_local=nothing)
    #=
    w1, w2, and w3 are weight factors for left moving particle f1, resting particle f0 (middle), and right moving particle, respectively  
    =#
    Q_eff = Q_local === nothing ? (@isdefined(Q) ? Q : 3) : Q_local
    D_eff = D_local === nothing ? (@isdefined(D) ? D : 1) : D_local

    @syms a, b, c, d

    if Q_eff == 3 && D_eff == 1
        @syms e1, e2, e3, w1, w2, w3
        w = [w1, w2, w3]
        w_value = [1. / 6, 2. / 3, 1. / 6]

        e = [e1, e2, e3]
        e_value = [-1.0, 0.0, 1.0]
    elseif Q_eff == 9 && D_eff == 2
        @syms w0 w1 w2 w3 w4 w5 w6 w7 w8
        @syms ex0 ey0 ex1 ey1 ex2 ey2 ex3 ey3 ex4 ey4 ex5 ey5 ex6 ey6 ex7 ey7 ex8 ey8

        w = [w0, w1, w2, w3, w4, w5, w6, w7, w8]
        w_value = [4. / 9, 1. / 9, 1. / 9, 1. / 9, 1. / 9, 1. / 36, 1. / 36, 1. / 36, 1. / 36]

        e = [
            [ex0, ey0],
            [ex1, ey1],
            [ex2, ey2],
            [ex3, ey3],
            [ex4, ey4],
            [ex5, ey5],
            [ex6, ey6],
            [ex7, ey7],
            [ex8, ey8],
        ]
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
        error("Unsupported lattice for lbm_const_sym: Q = $Q_eff, D = $D_eff")
    end

    a_value, b_value, c_value, d_value = 1., 3., 9. / 2., - 3. / 2.

    return w, e, w_value, e_value, a, b, c, d, a_value, b_value, c_value, d_value
end
