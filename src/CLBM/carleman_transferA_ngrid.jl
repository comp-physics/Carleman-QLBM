QCFD_SRC = ENV["QCFD_SRC"]  
QCFD_HOME = ENV["QCFD_HOME"]  

include(QCFD_SRC * "CLBM/carleman_transferA.jl")


function cal_delta_alpha(which_alpha, ngrid)
    delta_alpha = zeros(ngrid)
    delta_alpha[which_alpha] = 1.0
    return delta_alpha
end

function coeff_LBM_Fi_xalpha(Fj, which_alpha, ngrid)
    # F^(k) is a Q × Q^k matrix of constant coefficients of the polynomial combinations of f_m
    delta_alpha = cal_delta_alpha(which_alpha, ngrid)
    F_xalpha = kron(Fj, delta_alpha)
    return F_xalpha
end

function coeff_LBM_Fi_ngrid(Q, j, f, omega, tau_value, ngrid)
    len_F_xalpha = ngrid * Q
    Fj_ngrid = zeros(ngrid * len_F_xalpha, Q^j)
    # F^(k) (xalpha) a nQ × Q^k matrix of constant coefficients of the polynomial combinations of f_m
    # F^(k) (x) is a n^2 Q × Q^k matrix of constant coefficients of the polynomial combinations of f_m
    Fj = F_carlemanOrder_collision(Q, j, f, omega, tau_value)
    for i = 1 : ngrid
        F_xalpha = coeff_LBM_Fi_xalpha(Fj, i, ngrid)
        ind_s = (i - 1) * len_F_xalpha + 1
        ind_e = i * len_F_xalpha
        Fj_ngrid[ind_s:ind_e,:] = F_xalpha
    end
    return Fj, Fj_ngrid
end

function get_coeff_LBM_Fi_ngrid(poly_order, Q, f, omega, tau_value, ngrid) 
    upper_order = 3
    if poly_order != upper_order
        error("poly_order must be equal to upper_order")
    else    
        F1, F1_ngrid = coeff_LBM_Fi_ngrid(Q, upper_order - 2, f, omega, tau_value, ngrid)
        F2, F2_ngrid = coeff_LBM_Fi_ngrid(Q, upper_order - 1, f, omega, tau_value, ngrid)
        F3, F3_ngrid = coeff_LBM_Fi_ngrid(Q, upper_order, f, omega, tau_value, ngrid)
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
    phi = []
    for i = 1:ngrid
        phi_alpha = coeff_LBM_Fi_xalpha(f, i, ngrid)
        phi = append!(phi, phi_alpha)
    end
    return phi
end

function get_S_Fj(S, ngrid)
    #=
    Since Fj is broadcasted to n-point, i.e., Fj(x) = Fj(x_alpha) * delta_alpha, we need to have a S * delta_alpha version with j = 1 to match the dimension n^2Q * n^2Q of phi.  
    =#
    S_Fj = coeff_LBM_Fi_xalpha(S, 1, ngrid)
    for i = 2:ngrid
        S_alpha = coeff_LBM_Fi_xalpha(S, i, ngrid)
        S_Fj = hcat(S_Fj, S_alpha)
    end
    return S_Fj
end
