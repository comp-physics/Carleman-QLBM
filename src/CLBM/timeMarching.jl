using SparseArrays

function timeMarching_collision(omega, f, f_ini, tau_value, e_value, dt,  n_time, l_plot)
    omega_sub = LBM_const_subs(omega, tau_value)
  #  LB = lambdify(omega_sub .+ f, f)
    LB = lambdify(omega_sub * dt .+ f, f)
    #
#    f_ini = f_ini_test()
    Q = length(omega)
    fT = zeros(Q, n_time)
    uT = zeros(n_time)

    fT[:, 1] = f_ini 
    _, uT[1] = lbm_u(e_value, f_ini) 
    #println("fT = ", fT)

    for nt = 2:n_time
#        fT_temp = LB(fT[1, nt-1], fT[2, nt-1], fT[3, nt-1]) + F0_random_forcing(Q, force_factor, w_value, e_value) * dt
        fT_temp = LB(fT[1, nt-1], fT[2, nt-1], fT[3, nt-1]) + F0 * dt
        fT[:, nt] = fT_temp 
        _, uT[nt] = lbm_u(e_value, fT_temp) 
    end
    #
    if l_plot
        fm_plot(fT, n_time, ".", "k", "")
    end
    #
    return fT, uT
end

# Sparse version of carleman_C function
function carleman_C_sparse(Q, truncation_order, poly_order, f, omega, tau_value, force_factor, w_value, e_value)
    ncol_zero_ini = 0 # do NOT change this.
    C_dim = carleman_C_dim(Q, truncation_order, ngrid)
    
    # Initialize arrays for sparse matrix construction
    I_indices = Int[]
    J_indices = Int[]
    values = Float64[]
    
    #--b(t) term---
    F0 = F0_random_forcing(Q, force_factor, w_value, e_value)
    bt = spzeros(C_dim)
    bt[1:Q] = F0
    
    # Build sparse matrix by collecting non-zero blocks
    for ind_row = 1:truncation_order
        for ind_col = 1:truncation_order
            # Only fill blocks that satisfy the sparsity condition
            if ind_col >= ind_row - 1 && ind_col <= ind_row + poly_order - 1
                ind_row_C, ind_col_C = carleman_C_block_dim(Q, ind_row, ind_col, ncol_zero_ini)
                
                # Get the block matrix
                A_block = carleman_transferA(ind_row, ind_col, Q, f, omega, tau_value, force_factor, w_value, e_value, F0, ngrid)
                
                # Convert block to sparse and extract indices/values
                A_sparse = sparse(A_block)
                block_I, block_J, block_vals = findnz(A_sparse)
                
                # Adjust indices to global matrix coordinates
                global_I = block_I .+ (first(ind_row_C) - 1)
                global_J = block_J .+ (first(ind_col_C) - 1)
                
                # Append to global arrays
                append!(I_indices, global_I)
                append!(J_indices, global_J)
                append!(values, block_vals)
            end
        end
    end
    
    # Create sparse matrix
    C_sparse = sparse(I_indices, J_indices, values, C_dim, C_dim)
    
    return C_sparse, bt, F0
end

#
function timeMarching_collision_CLBM(omega, f, tau_value, Q, C, truncation_order, e_value, dt, f_ini, n_time, l_plot)
    V0 = carleman_V(f_ini, truncation_order)
    V0 = Float64.(V0)
#    CL = C .* dt .+ Matrix{Float64}(I, size(C)) # explicit Euler scheme


    VT = zeros(size(C)[1], n_time)
    VT[:, 1] = V0

    VT_f = zeros(Q, n_time)
    VT_f[:, 1] = VT[1:Q, 1] 

    uT = zeros(n_time)
    _, uT[1] = lbm_u(e_value, VT_f[:, 1]) 
#    println("VT_f=", VT_f)
    #---LBM---
    omega_sub = LBM_const_subs(omega, tau_value)
    LB = lambdify(omega_sub * dt .+ f, f)
    #
    fT = zeros(Q, n_time)
    fT[:, 1] = f_ini 
    #---
    #---time marching---
    for nt = 2:n_time
        C, bt, F0 = carleman_C(Q, truncation_order, poly_order, f, omega, tau_value, force_factor, w_value, e_value)
#        CL = C .* dt .+ Matrix{Float64}(I, size(C))
#        VT[:, nt] = CL * VT[:, nt - 1]
        VT[:, nt] = (C * VT[:, nt - 1] + bt) .* dt .+ VT[:, nt - 1]
        _, uT[nt] = lbm_u(e_value, VT[1:Q, nt]) 
        #---LBM---
        fT_temp = LB(fT[1, nt-1], fT[2, nt-1], fT[3, nt-1]) + F0 * dt
        fT[:, nt] = fT_temp 
    end
    VT_f = VT[1:Q, :]
    #
    if l_plot
        fm_plot(VT_f, n_time, ".", "k", "")
    end
    #
    return VT_f, VT, uT, fT
end

# Sparse version of the time marching function
function timeMarching_collision_CLBM_sparse(omega, f, tau_value, Q, truncation_order, e_value, dt, f_ini, n_time, l_plot)
    V0 = carleman_V(f_ini, truncation_order)
    V0 = Float64.(V0)

    # Initialize with proper size based on V0 length
    VT = zeros(length(V0), n_time)
    VT[:, 1] = V0

    VT_f = zeros(Q, n_time)
    VT_f[:, 1] = VT[1:Q, 1] 

    uT = zeros(n_time)
    _, uT[1] = lbm_u(e_value, VT_f[:, 1]) 

    #---LBM---
    omega_sub = LBM_const_subs(omega, tau_value)
    LB = lambdify(omega_sub * dt .+ f, f)
    
    fT = zeros(Q, n_time)
    fT[:, 1] = f_ini 
    
    #---time marching---
    for nt = 2:n_time
        # Use sparse matrix version
        C_sparse, bt, F0 = carleman_C_sparse(Q, truncation_order, poly_order, f, omega, tau_value, force_factor, w_value, e_value)
        
        # Sparse matrix-vector multiplication
        VT[:, nt] = (C_sparse * VT[:, nt - 1] + bt) .* dt .+ VT[:, nt - 1]
        _, uT[nt] = lbm_u(e_value, VT[1:Q, nt]) 
        
        #---LBM---
        fT_temp = LB(fT[1, nt-1], fT[2, nt-1], fT[3, nt-1]) + F0 * dt
        fT[:, nt] = fT_temp 
    end
    VT_f = VT[1:Q, :]
    
    if l_plot
        fm_plot(VT_f, n_time, ".", "k", "")
    end
    
    return VT_f, VT, uT, fT
end
