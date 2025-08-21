QCFD_SRC = ENV["QCFD_SRC"]  
include(QCFD_SRC * "LBM/cal_feq.jl")

function f_ini_feq(u0, rho0)
    _, _, w_value, e_value, _, _, _, _, a_v, b_v, c_v, d_v = lbm_const_sym()
    vx = u0; vy = 0.; wm = w_value; vcx = e_value; vcy = [0. , 0. , 0.]
    f_ini = cal_feq(rho0, vx, vy, wm, vcx, vcy, a_v, b_v, c_v, d_v)
    println("sum f_m =", sum(f_ini))
    return f_ini
end

function f_ini_test(uu)
    _, _, w_value, _, _, _, _, _, a_v, b_v, c_v, d_v = lbm_const_sym()
    f_ini = w_value .+ [-uu/2, 0, uu/2]
    return f_ini 
end
