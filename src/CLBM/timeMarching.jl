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

# Optimized sparse version of carleman_C function
function carleman_C_sparse(Q, truncation_order, poly_order, f, omega, tau_value, force_factor, w_value, e_value)
    ncol_zero_ini = 0 # do NOT change this.
    C_dim = carleman_C_dim(Q, truncation_order, ngrid)
    
    # OPTIMIZATION: Estimate non-zero elements for better memory pre-allocation
    estimated_nnz = estimate_carleman_nnz(Q, truncation_order, poly_order, ngrid)
    
    # Pre-allocate arrays with estimated size to reduce reallocations
    I_indices = Vector{Int}()
    J_indices = Vector{Int}()
    values = Vector{Float64}()
    sizehint!(I_indices, estimated_nnz)
    sizehint!(J_indices, estimated_nnz)
    sizehint!(values, estimated_nnz)
    
    #--b(t) term---
    F0 = F0_random_forcing(Q, force_factor, w_value, e_value)
    bt = spzeros(C_dim)
    bt[1:Q] = F0
    
    # MEMORY-EFFICIENT: Build sparse matrix incrementally without storing all positions
    # For large matrices, sacrifice perfect overlap handling for memory efficiency
    
    # Estimate memory requirement first
    dense_memory_gb = C_dim^2 * 8 / 1024^3
    
    if dense_memory_gb > 5.0
        # LARGE MATRIX: Use accumulation method (may not handle overlaps perfectly)
        println("⚠️  Large matrix detected ($(round(dense_memory_gb, digits=1)) GB dense)")
        println("   Using memory-efficient sparse assembly...")
        
        I_indices = Vector{Int}()
        J_indices = Vector{Int}()
        values = Vector{Float64}()
        
        for ind_row = 1:truncation_order
            for ind_col = 1:truncation_order
                if ind_col >= ind_row - 1 && ind_col <= ind_row + poly_order - 1
                    ind_row_C, ind_col_C = carleman_C_block_dim(Q, ind_row, ind_col, ncol_zero_ini)
                    A_block_sparse = carleman_transferA_sparse(ind_row, ind_col, Q, f, omega, tau_value, force_factor, w_value, e_value, F0, ngrid)
                    
                    # Extract and add non-zeros directly
                    block_I, block_J, block_vals = findnz(A_block_sparse)
                    row_offset = first(ind_row_C) - 1
                    col_offset = first(ind_col_C) - 1
                    
                    for k = 1:length(block_I)
                        push!(I_indices, block_I[k] + row_offset)
                        push!(J_indices, block_J[k] + col_offset)
                        push!(values, block_vals[k])
                    end
                end
            end
        end
        
        # Build sparse matrix (Julia handles duplicate indices by summing - approximate)
        if !isempty(I_indices)
            C_sparse = sparse(I_indices, J_indices, values, C_dim, C_dim)
            println("   Matrix construction completed: $(nnz(C_sparse)) non-zeros")
        else
            C_sparse = spzeros(C_dim, C_dim)
        end
        
    else
        # SMALL MATRIX: Use exact method with dictionary for perfect correctness
        final_values = Dict{Tuple{Int,Int}, Float64}()
        
        for ind_row = 1:truncation_order
            for ind_col = 1:truncation_order
                if ind_col >= ind_row - 1 && ind_col <= ind_row + poly_order - 1
                    ind_row_C, ind_col_C = carleman_C_block_dim(Q, ind_row, ind_col, ncol_zero_ini)
                    A_block_sparse = carleman_transferA_sparse(ind_row, ind_col, Q, f, omega, tau_value, force_factor, w_value, e_value, F0, ngrid)
                    
                    # Convert to dense block for exact overwrite semantics
                    A_block_dense = Array(A_block_sparse)
                    
                    # Store all values (including zeros) for exact semantics
                    for i = 1:size(A_block_dense, 1), j = 1:size(A_block_dense, 2)
                        global_row = ind_row_C[i]
                        global_col = ind_col_C[j]
                        final_values[(global_row, global_col)] = A_block_dense[i, j]
                    end
                end
            end
        end
        
        # Build sparse matrix from non-zero final values
        I_final = Int[]
        J_final = Int[]
        vals_final = Float64[]
        
        for ((row, col), value) in final_values
            if abs(value) > 1e-15
                push!(I_final, row)
                push!(J_final, col)
                push!(vals_final, value)
            end
        end
        
        if !isempty(I_final)
            C_sparse = sparse(I_final, J_final, vals_final, C_dim, C_dim)
        else
            C_sparse = spzeros(C_dim, C_dim)
        end
    end
    
    return C_sparse, bt, F0
end

# OPTIMIZATION: Function to estimate non-zero elements for better memory allocation
function estimate_carleman_nnz(Q, truncation_order, poly_order, ngrid)
    """
    Estimate the number of non-zero elements in the Carleman matrix for better memory pre-allocation.
    This reduces the number of reallocations during sparse matrix assembly.
    """
    total_blocks = 0
    
    # Count the number of active blocks based on sparsity pattern
    for ind_row = 1:truncation_order
        for ind_col = 1:truncation_order
            if ind_col >= ind_row - 1 && ind_col <= ind_row + poly_order - 1
                total_blocks += 1
            end
        end
    end
    
    # Estimate average non-zeros per block
    # For Kronecker products of sparse matrices, density typically scales as:
    # - Q^i matrices have roughly Q*i non-zeros for collision operators
    # - Identity kronecker products maintain sparsity well
    if ngrid == 1
        # For ngrid=1, blocks are roughly Q^i × Q^j in size
        avg_block_nnz = Q * truncation_order * 2  # Conservative estimate
    else
        # For ngrid>1, blocks grow as (Q*ngrid^2)^i but sparsity increases
        avg_block_nnz = Q * ngrid * truncation_order  # More conservative for larger grids
    end
    
    estimated_nnz = total_blocks * avg_block_nnz
    
    # Add 20% buffer to account for estimation uncertainty
    return Int(ceil(estimated_nnz * 1.2))
end

#
function timeMarching_collision_CLBM(omega, f, tau_value, Q, C, truncation_order, e_value, dt, f_ini, n_time, l_plot)
    # replace f_ini with a well-developed snapshot of stochastic D1Q3 LBM: dt1_force0.00025_tau0.503_Nx16_NT50000010.h5.
    V0 = carleman_V(f_ini, truncation_order)
    V0 = Float64.(V0)

    # FIXED: Use the pre-computed matrix C passed as parameter (built outside this function)
    # Also need to compute bt and F0 once (not in the time loop)
    _, bt, F0 = carleman_C(Q, truncation_order, poly_order, f, omega, tau_value, force_factor, w_value, e_value)

    VT = zeros(size(C)[1], n_time)
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
        # FIXED: Use passed matrix C and computed bt, F0 - no reconstruction needed
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

    # FIXED: Build sparse matrix ONCE outside the time loop for efficiency
    C_sparse, bt, F0 = carleman_C_sparse(Q, truncation_order, poly_order, f, omega, tau_value, force_factor, w_value, e_value)

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
        # FIXED: Use pre-computed sparse matrix - no reconstruction needed
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

# ========================================================================
# SPARSE KRONECKER PRODUCT FUNCTIONS FOR MEMORY-EFFICIENT OPERATIONS
# ========================================================================

function Kron_kth_sparse(ff, k)
    if k == 1
        return issparse(ff) ? ff : sparse(ff)
    else
        # Ensure input is sparse from the start
        ff_sparse = issparse(ff) ? ff : sparse(ff)
        fk = ff_sparse
        for i = 1:k-1
            # FIXED: Match the order in original Kron_kth: kron(ff, fk)
            fk = kron(ff_sparse, fk)
        end
        return fk
    end
end

function Kron_kth_identity_sparse(Fj, i, rth, Q)
    if !@isdefined(cached_sparse_identity) || size(cached_sparse_identity, 1) != Q
        global cached_sparse_identity = sparse(1.0I, Q, Q)
    end
    identity_sparse = cached_sparse_identity
    
    if rth > i
        error("rth must be smaller than i")
    end
    
    if i == 1
        return issparse(Fj) ? Fj : sparse(Fj)
    else
        # Ensure Fj is sparse once at the beginning
        Fj_sparse = issparse(Fj) ? Fj : sparse(Fj)
        
        if rth == 1
            imatrix_right = Kron_kth_sparse(identity_sparse, i - rth)
            return kron(Fj_sparse, imatrix_right)
        elseif rth == i
            imatrix_left = Kron_kth_sparse(identity_sparse, rth - 1)
            return kron(imatrix_left, Fj_sparse)
        else
            imatrix_left = Kron_kth_sparse(identity_sparse, rth - 1)
            imatrix_right = Kron_kth_sparse(identity_sparse, i - rth)
            A_sub = kron(imatrix_left, Fj_sparse)
            return kron(A_sub, imatrix_right)
        end
    end
end

function sum_Kron_kth_identity_sparse(Fj, i, Q)
    """Sparse version of sum_Kron_kth_identity"""
    A_ij = Kron_kth_identity_sparse(Fj, i, 1, Q)
    for rth = 2:i
        A_ij = A_ij + Kron_kth_identity_sparse(Fj, i, rth, Q)
    end
    return A_ij
end

function transferA_ngrid_sparse(i, j, Q, ngrid)
    """Sparse version of transferA_ngrid - avoids creating large dense matrices"""
    # Get the appropriate F matrix
    if j == 1
        Fj_ngrid = F1_ngrid
    elseif j == 2
        Fj_ngrid = F2_ngrid
    elseif j == 3
        Fj_ngrid = F3_ngrid
    else
        error("j of F^{j} must be 1, 2, 3, ..., poly_order")
    end
    
    # Use sparse version of sum_Kron_kth_identity
    A_ij = sum_Kron_kth_identity_sparse(Fj_ngrid, i, Q * ngrid)
    
    return A_ij
end

function carleman_transferA_sparse(ind_row, ind_col, Q, f, omega, tau_value, force_factor, w_value, e_value, F0, ngrid)
    """Sparse version of carleman_transferA"""
    if ind_row <= ind_col 
        i = ind_row
        if ngrid > 1
            j = ind_col
        else
            j = Int(ind_col - (i - 1))
        end
        A = transferA_ngrid_sparse(i, j, Q, ngrid)
    else
       # The A_{i+j-1}^i with i >= 1 and j = 0
        i = ind_row 
        j = i - 1 
        if ngrid > 1
            A = transferA_ngrid_sparse(i, j, Q, ngrid)
            A = spzeros(size(A)...)  # Return sparse zero matrix
        else
            # Use sparse version for F0 term
            A = sum_Kron_kth_identity_sparse(F0, i, Q * ngrid)
        end
    end
    
    return A
end


