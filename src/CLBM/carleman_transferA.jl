QCFD_HOME = ENV["QCFD_HOME"]  
QCFD_SRC = ENV["QCFD_SRC"]  
include(QCFD_HOME * "/julia_lib/matrix_kit.jl")
include(QCFD_SRC * "CLBM/LBM_const_subs.jl")
#include(QCFD_SRC * "LBM/julia/forcing.jl")

function F_carlemanOrder_Q_collision(iQ, order, f, omega, tau_value)
    #=
    1. This function is to calculate the coefficient of the collision matrix at each row
    2. somehow the error "UndefVarError: `coeffs` not defined" occurred when if condition is used here
    =#
#    if l_sympy
        vars, coll_coeff = coeffs(omega[iQ], f)
#    else
#        coeffs, _ = polynomial_coeffs(omega[iQ], f)
#        vars = collect(keys(coeffs))  
#        coll_coeff = collect(values(coeffs)) 
#    end

    if l_sympy
        ind = findall(sum([degree(i, j) for i in vars, j in f], dims=(2,3)) .== order)
    else
        ind = findall([Symbolics.degree(i) == order for i in vars])
    end
    vars_F2 = vars[ind]
    coll_coeff_F2 = coll_coeff[ind] #XY: to substitute vars
    coll_coeff_F2 = LBM_const_subs(coll_coeff_F2, tau_value)

    if l_sympy
       # f_kron = kron(f, f)
        f_kron = Kron_kth(f, order)
    else
        f_kron = kronecker(f, f)
    end

    coeff_kron = Array{Float64}(undef, length(f_kron))
    coeff_kron[:] .= 0.

    if l_sympy
        for i = 1:length(vars_F2)
            local ind_temp
            ind_temp = findall(f_kron .== vars_F2[i])[1]
           coeff_kron[ind_temp] = Float64.(coll_coeff_F2[i])
        end
    else
        for i = 1:length(vars_F2)
            local ind_temp
            ind_temp = findall(isequal.(f_kron, vars_F2[i]))[1]
            coeff_kron[ind_temp] = coll_coeff_F2[i]
        end
    end
    #
    return f_kron, coeff_kron
end
#
function F_carlemanOrder_collision(Q, order, f, omega, tau_value)
    F_k = Array{Float64}(undef, (Q, Int(Q^order)))
    for i = 1:Q
        f_kron, coeff_kron = F_carlemanOrder_Q_collision(i, order, f, omega, tau_value)
        F_k[i, :] = coeff_kron   
    end
    #
    return F_k
end

function F0_random_forcing(Q, force_factor, w_value, e_value)
    LX = 1
    LY = 1
    j = 1
    k = 1
    fx_domain, fy_domain = forcing_random(LX, LY, 1, force_factor)
    F0 = forcing_lbm(fx_domain[j, k], fy_domain[j, k], w_value, e_value)
    return F0
end

function Kron_kth_identity(Fj, i, rth, Q)
    #
    identity_matrix = Matrix(1.0I, Q, Q)
    #
    if rth > i
        exit("rth must be smaller than i")
    end
    #
    if i == 1
        A_sub = Fj
    else
        if rth == 1
            imatrix_right = Kron_kth(identity_matrix, i - rth) 
            A_sub = kron(Fj, imatrix_right)
        elseif rth == i
            imatrix_left = Kron_kth(identity_matrix, rth - 1)
            A_sub = kron(imatrix_left, Fj)
        else
            imatrix_left = Kron_kth(identity_matrix, rth - 1)
            imatrix_right = Kron_kth(identity_matrix, i - rth )
            A_sub = kron(imatrix_left, Fj)
            A_sub = kron(A_sub, imatrix_right)
        end
    end
    #
    return A_sub
    #
end

function sum_Kron_kth_identity(Fj, i, Q)
    A_ij = Kron_kth_identity(Fj, i, 1, Q)
    for rth = 2:i
        A_ij = A_ij + Kron_kth_identity(Fj, i, rth, Q) 
    end
    return A_ij
end

function transferA(i, j, Q, f, omega, tau_value)
    # the transfer matrix A_i^{i + j - 1} has a dimension of Q^i by Q^{i + j -1} for given index i and j 
    Fj = F_carlemanOrder_collision(Q, j, f, omega, tau_value)
    A_ij = sum_Kron_kth_identity(Fj, i, Q)
    return A_ij
end

#function transferA_F0(i, j, Q, force_factor, w_value, e_value)
function transferA_F0(i, Q, force_factor, w_value, e_value, F0, ngrid)
    # the transfer matrix A_i^{i + j - 1} has a dimension of Q^i by Q^{i + j -1} for given index i and j 
#    Fj = F0_random_forcing(Q, force_factor, w_value, e_value)
#    A_ij = sum_Kron_kth_identity(Fj, i, Q)
    A_ij = sum_Kron_kth_identity(F0, i, Q * ngrid)
    return A_ij
end

function carleman_transferA(ind_row, ind_col, Q, f, omega, tau_value, force_factor, w_value, e_value, F0, ngrid)
    if ind_row <= ind_col 
        i = ind_row
        if ngrid > 1
            j = ind_col
        else
            j = Int(ind_col - (i - 1))
        end
        A = transferA_ngrid(i, j, Q, ngrid)
    else
       # exit("ind_row must be <= ind_col")
       # The A_{i+j-1}^i with i >= 1 and j = 0
        i = ind_row 
        j = i - 1 
        if ngrid > 1
            A = transferA_ngrid(i, j, Q, ngrid)
            A[:] .= 0.0
        else
            A = transferA_F0(i, Q, force_factor, w_value, e_value, F0, ngrid)
        end
    end
    #
    return A
end

function carleman_transferA_S(ind_row, ind_col, Q, ngrid, S_Fj)
    if ind_row == ind_col 
        i = ind_row
        A = transferA_S(i, Q, ngrid, S_Fj)
    end
    #
    return A
end


function carleman_C_dim(Q, truncation_order, ngrid)
    C_dim = 0
    for i = 1:truncation_order
       # C_dim = C_dim + (Q * ngrid) ^ i
       # C_dim = C_dim + ngrid * (Q * ngrid) ^ i
        C_dim = C_dim + (Q * ngrid^2) ^ i
    end
    return C_dim
end

function A_row_dim(Q, truncation_order, ngrid)
    C_dim = 0
    for i = 1:truncation_order
       # C_dim = C_dim + (Q * ngrid) ^ i
        C_dim = C_dim + ngrid * (Q * ngrid) ^ i
    end
    return C_dim
end


function A_col_dim(Q, i, j, ngrid)

    C_dim = ngrid ^ (i-1) * Q ^ (i + j -1)
    return C_dim
end


function carleman_C_block_dim(Q, ind_row, ind_col, ncol_zero_ini)
    if ind_row == 1 && ind_col == 1
#        ind_row_C = Int((Q * ngrid) ^ (ind_row - 1)):Int(carleman_C_dim(Q, ind_row, ngrid)) 
#        ind_col_C = Int(Q ^ (ind_col - 1)):Int(carleman_C_dim(Q, ind_col, ngrid)) 
        ind_row_C = 1:Int(carleman_C_dim(Q, ind_row, ngrid)) 
#        ind_col_C = 1:Int(carleman_C_dim(Q, ind_col, ngrid)) 
        ind_col_C = 1:Int(A_col_dim(Q, ind_row, ind_col, ngrid)) 
    else
        if ngrid > 1
            ind_row_C = Int(A_row_dim(Q, ind_row - 1, ngrid) + 1):Int(A_row_dim(Q, ind_row, ngrid)) 
            ind_col_C = Int(A_col_dim(Q, ind_row, ind_col - 1, ngrid) + 1 - ncol_zero_ini):Int(A_col_dim(Q, ind_row, ind_col - 1, ngrid) + A_col_dim(Q, ind_row, ind_col, ngrid) - ncol_zero_ini) 
        else
            ind_row_C = Int(carleman_C_dim(Q, ind_row - 1, ngrid) + 1):Int(carleman_C_dim(Q, ind_row, ngrid)) 
            ind_col_C = Int(carleman_C_dim(Q, ind_col - 1, ngrid) + 1 - ncol_zero_ini):Int(carleman_C_dim(Q, ind_col, ngrid) - ncol_zero_ini) 
        end
    end

    return ind_row_C, ind_col_C 
end

function carleman_C(Q, truncation_order, poly_order, f, omega, tau_value, force_factor, w_value, e_value)
    ncol_zero_ini = 0 # do NOT change this.
    C_dim = carleman_C_dim(Q, truncation_order, ngrid)
    C = Array{Float64}(undef, (C_dim, C_dim))
    C[:] .= 0.
    #--b(t) term---
    #
    F0 = F0_random_forcing(Q, force_factor, w_value, e_value)
    #
    bt = zeros(C_dim)
   # bt[1:Q] = F0_random_forcing(Q, force_factor, w_value, e_value)
    bt[1:Q] = F0
    #
    for ind_row = 1:truncation_order
        for ind_col = 1:truncation_order
            # if ind_col >= ind_row && ind_col <= ind_row + poly_order - 1 #A_{i+j-1}^i with i >= 1 and j >= 1 
            if ind_col >= ind_row - 1 && ind_col <= ind_row + poly_order - 1 #A_{i+j-1}^i with i >= 1 and j = 0
                ind_row_C, ind_col_C = carleman_C_block_dim(Q, ind_row, ind_col, ncol_zero_ini)
                C[ind_row_C, ind_col_C] = carleman_transferA(ind_row, ind_col, Q, f, omega, tau_value, force_factor, w_value, e_value, F0, ngrid)
            end
        end
    end
    #
    return C, bt, F0
end

function carleman_S(Q, truncation_order, poly_order, ngrid, S_Fj)
    ncol_zero_ini = 0 # do NOT change this.
    C_dim = carleman_C_dim(Q, truncation_order, ngrid)
    S = Array{Float64}(undef, (C_dim, C_dim))
    S[:] .= 0.
    #
    for ind_row = 1:truncation_order
        for ind_col = 1:truncation_order
            if ind_col >= ind_row - 1 && ind_col <= ind_row + poly_order - 1
                if ind_col == ind_row
                    ind_row_C, _ = carleman_C_block_dim(Q, ind_row, ind_col, ncol_zero_ini)
                    S[ind_row_C, ind_row_C] = carleman_transferA_S(ind_row, ind_col, Q, ngrid, S_Fj)
                end
            end
        end
    end
    #
    return S
end


function carleman_V(f, truncation_order)
    V = []
    for i = 1:truncation_order
        V = append!(V, Kron_kth(f, i))
    end
    #
    return V
end
