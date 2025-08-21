using SparseArrays
## Function to create the gradient matrix with periodic boundary conditions
#function gradient_matrix(n::Int, h::Float64)
#    # Initialize an n×n matrix of zeros
#    D = zeros(Float64, n, n)
#    
#    # Fill the matrix according to central difference with periodic BC
#    for i in 1:n
#        if i == 1
#            # First row wraps around to the last element
#            D[i, n] = -1
#            D[i, i+1] = 1
#        elseif i == n
#            # Last row wraps around to the first element
#            D[i, i-1] = -1
#            D[i, 1] = 1
#        else
#            # Central difference for the internal elements
#            D[i, i-1] = -1
#            D[i, i+1] = 1
#        end
#    end
#    
#    # Scale the matrix by 1 / (2h)
#    D *= 1 / (2 * h)
#    
#    return D
#end
#
## Example usage:
#h = 1.0       # Grid spacing
#D = gradient_matrix(ngrid, h)
#
#println("Gradient Matrix D with Periodic Boundary Conditions:")
#println(D)
#
#function get_em_ngrid(e_value, m, ngrid)
#    em_ngrid = zeros(ngrid)
#    em_ngrid[:] .= e_value[m]
#    return em_ngrid
#end
#
#e1_ngrid = get_em_ngrid(e_value, 1, ngrid)
#S1 = e1_ngrid' * D 
#
#
#function get_S_ngrid(e_value, Q,  ngrid, D)
#    S_ngrid = zeros(Q, ngrid)
#    for m = 1:Q
#        em_ngrid = get_em_ngrid(e_value, m, ngrid)
#        S_ngrid[m, :] = em_ngrid' * D
#    end
#    return S_ngrid
#end
#
#S_ngrid = get_S_ngrid(e_value, Q,  ngrid, D)

# w, e, w_value, e_value, a, b, c, d, a_value, b_value, c_value, d_value = lbm_const_sym()

# Function to create the 1D LBM streaming matrix
function streaming_matrix(n::Int)
    # Define the base blocks for f2 (no movement), f3 (right movement), f1 (left movement)
    S2 = Matrix{Float64}(1.0I, n, n)  # Identity for f2 (no movement)

    # Right shift block with periodic boundary conditions
    S3 = circshift(S2, (0, -1))

    # Left shift block with periodic boundary conditions
    S1 = circshift(S2, (0, 1))
    #= To match the definition of phi(x), return S' =  
    [f1(x1), f2(x1), ..., fQ(x1)];
    [f1(x2), f2(x2), ..., fQ(x2)];
    .
    [f1(xn), f2(xn), ..., fQ(xn)];
    =# 
    # Combine the blocks into the full streaming matrix
    S = 
         [zeros(n, n)  zeros(n, n) S1;
         S2  zeros(n, n)  zeros(n, n);
         zeros(n, n)  S3  zeros(n, n)]
         
    return S'
end



function sparse_central_diff(n, h)
    """
    Create sparse tri-diagonal matrix for large systems
    """
    
    # Row and column indices
    I = Int[]
    J = Int[]
    V = Float64[]
    
    for i in 1:n
        # Main diagonal
        push!(I, i)
        push!(J, i)
        push!(V, -2.0/h^2)
        
        # Super diagonal
        if i < n
            push!(I, i)
            push!(J, i+1)
            push!(V, 1.0/h^2)
        end
        
        # Sub diagonal
        if i > 1
            push!(I, i)
            push!(J, i-1)
            push!(V, 1.0/h^2)
        end
    end
    
    return sparse(I, J, V, n, n)
end

# Large system example
n = 4
#h = 1.0/(n+1)
h = 1.0
A_sparse = sparse_central_diff(n, h)
println("Sparse matrix size: $(size(A_sparse))")
println("Non-zero elements: $(nnz(A_sparse))")

function streaming_operator_D1Q3(nx, hx)
    """
    Construct streaming operators for D1Q3 lattice
    
    Args:
        nx: number of grid points in x-direction
        hx: grid spacing in x-direction
    
    Returns:
        Array of sparse matrices [S0, S1, S2] for each velocity
    """
    
    # D1Q3 velocity vectors
    e = [
        [0],    # e0: rest particle
        [1],    # e1: moving right (+x)
        [-1]    # e2: moving left (-x)
    ]
    
    streaming_ops = []
    
    for i in 1:3
        ex = e[i][1]
        
        # Create sparse matrix for ∂/∂x
        I_idx = Int[]
        J_idx = Int[]
        V_vals = Float64[]
        
        for j in 1:nx
            if ex != 0  # Non-zero velocity component
                # Central difference: ∂f/∂x ≈ (f[j+1] - f[j-1])/(2*hx)
                if j > 1  # Left neighbor exists
                    push!(I_idx, j)
                    push!(J_idx, j-1)
                    push!(V_vals, -ex/(2*hx))
                end
                
                if j < nx  # Right neighbor exists
                    push!(I_idx, j)
                    push!(J_idx, j+1)
                    push!(V_vals, ex/(2*hx))
                end
            else
                # Rest particle: no streaming
                # Add zero diagonal (optional, for completeness)
                push!(I_idx, j)
                push!(J_idx, j)
                push!(V_vals, 0.0)
            end
        end
        
        S = sparse(I_idx, J_idx, V_vals, nx, nx)
        push!(streaming_ops, S)
    end
    S_combined = hcat(streaming_ops...)
    # Combine all streaming operators into a single matrix

    return streaming_ops, e, S_combined
end

# Example usage D1Q3
nx = 5
hx = 1
S_D1Q3, e_D1Q3, streaming_matrix_out = streaming_operator_D1Q3(nx, hx)

println("D1Q3 Streaming Operators:")
println("=" ^ 40)
for i in 1:3
    println("Velocity e[$i] = $(e_D1Q3[i])")
    println("Streaming operator S[$i]:")
    println(Matrix(S_D1Q3[i]))
    println()
end

function streaming_operator_D2Q9(nx, ny, hx, hy)
    """
    Construct streaming operators for D2Q9 lattice
    
    Args:
        nx, ny: number of grid points in x, y directions
        hx, hy: grid spacing in x, y directions
    
    Returns:
        Array of sparse matrices for each velocity
    """
    
    # D2Q9 velocity vectors
    e = [
        [0, 0],   # e0: rest
        [1, 0],   # e1: +x
        [0, 1],   # e2: +y  
        [-1, 0],  # e3: -x
        [0, -1],  # e4: -y
        [1, 1],   # e5: +x, +y
        [-1, 1],  # e6: -x, +y
        [-1, -1], # e7: -x, -y
        [1, -1]   # e8: +x, -y
    ]
    
    n_total = nx * ny
    streaming_ops = []
    
    # Helper function to convert 2D indices to linear index
    idx_2d_to_1d(i, j) = (j-1) * nx + i
    
    for vel_idx in 1:9
        ex, ey = e[vel_idx]
        
        I_idx = Int[]
        J_idx = Int[]
        V_vals = Float64[]
        
        for j in 1:ny
            for i in 1:nx
                linear_idx = idx_2d_to_1d(i, j)
                
                # x-direction streaming: ex * ∂/∂x
                if ex != 0
                    if i > 1  # Left neighbor
                        neighbor_idx = idx_2d_to_1d(i-1, j)
                        push!(I_idx, linear_idx)
                        push!(J_idx, neighbor_idx)
                        push!(V_vals, -ex / (2*hx))
                    end
                    
                    if i < nx  # Right neighbor
                        neighbor_idx = idx_2d_to_1d(i+1, j)
                        push!(I_idx, linear_idx)
                        push!(J_idx, neighbor_idx)
                        push!(V_vals, ex / (2*hx))
                    end
                end
                
                # y-direction streaming: ey * ∂/∂y
                if ey != 0
                    if j > 1  # Bottom neighbor
                        neighbor_idx = idx_2d_to_1d(i, j-1)
                        push!(I_idx, linear_idx)
                        push!(J_idx, neighbor_idx)
                        push!(V_vals, -ey / (2*hy))
                    end
                    
                    if j < ny  # Top neighbor
                        neighbor_idx = idx_2d_to_1d(i, j+1)
                        push!(I_idx, linear_idx)
                        push!(J_idx, neighbor_idx)
                        push!(V_vals, ey / (2*hy))
                    end
                end
                
                # Rest particle case
                if ex == 0 && ey == 0
                    push!(I_idx, linear_idx)
                    push!(J_idx, linear_idx)
                    push!(V_vals, 0.0)
                end
            end
        end
        
        S = sparse(I_idx, J_idx, V_vals, n_total, n_total)
        push!(streaming_ops, S)
    end
    S_combined = hcat(streaming_ops...)
    return streaming_ops, e, S_combined
end

function streaming_operator_D3Q27(nx, ny, nz, hx, hy, hz)
    """
    Construct streaming operators for D3Q27 lattice
    
    Args:
        nx, ny, nz: number of grid points in x, y, z directions
        hx, hy, hz: grid spacing in x, y, z directions
    
    Returns:
        Array of sparse matrices for each velocity
    """
    
    # D3Q27 velocity vectors
    e = []
    
    # Generate all combinations of {-1, 0, 1}³
    for ez in [-1, 0, 1]
        for ey in [-1, 0, 1]
            for ex in [-1, 0, 1]
                push!(e, [ex, ey, ez])
            end
        end
    end
    
    n_total = nx * ny * nz
    streaming_ops = []
    
    # Helper function to convert 3D indices to linear index
    idx_3d_to_1d(i, j, k) = (k-1) * nx * ny + (j-1) * nx + i
    
    for vel_idx in 1:27
        ex, ey, ez = e[vel_idx]
        
        I_idx = Int[]
        J_idx = Int[]
        V_vals = Float64[]
        
        for k in 1:nz
            for j in 1:ny
                for i in 1:nx
                    linear_idx = idx_3d_to_1d(i, j, k)
                    
                    # x-direction streaming: ex * ∂/∂x
                    if ex != 0
                        if i > 1  # Left neighbor
                            neighbor_idx = idx_3d_to_1d(i-1, j, k)
                            push!(I_idx, linear_idx)
                            push!(J_idx, neighbor_idx)
                            push!(V_vals, -ex / (2*hx))
                        end
                        
                        if i < nx  # Right neighbor
                            neighbor_idx = idx_3d_to_1d(i+1, j, k)
                            push!(I_idx, linear_idx)
                            push!(J_idx, neighbor_idx)
                            push!(V_vals, ex / (2*hx))
                        end
                    end
                    
                    # y-direction streaming: ey * ∂/∂y
                    if ey != 0
                        if j > 1  # Bottom neighbor
                            neighbor_idx = idx_3d_to_1d(i, j-1, k)
                            push!(I_idx, linear_idx)
                            push!(J_idx, neighbor_idx)
                            push!(V_vals, -ey / (2*hy))
                        end
                        
                        if j < ny  # Top neighbor
                            neighbor_idx = idx_3d_to_1d(i, j+1, k)
                            push!(I_idx, linear_idx)
                            push!(J_idx, neighbor_idx)
                            push!(V_vals, ey / (2*hy))
                        end
                    end
                    
                    # z-direction streaming: ez * ∂/∂z
                    if ez != 0
                        if k > 1  # Back neighbor
                            neighbor_idx = idx_3d_to_1d(i, j, k-1)
                            push!(I_idx, linear_idx)
                            push!(J_idx, neighbor_idx)
                            push!(V_vals, -ez / (2*hz))
                        end
                        
                        if k < nz  # Front neighbor
                            neighbor_idx = idx_3d_to_1d(i, j, k+1)
                            push!(I_idx, linear_idx)
                            push!(J_idx, neighbor_idx)
                            push!(V_vals, ez / (2*hz))
                        end
                    end
                    
                    # Rest particle case
                    if ex == 0 && ey == 0 && ez == 0
                        push!(I_idx, linear_idx)
                        push!(J_idx, linear_idx)
                        push!(V_vals, 0.0)
                    end
                end
            end
        end
        
        S = sparse(I_idx, J_idx, V_vals, n_total, n_total)
        push!(streaming_ops, S)
    end
    S_combined = hcat(streaming_ops...)
    return streaming_ops, e, S_combined
end

function streaming_operator_D1Q3_interleaved(nx, hx)
    """
    Create D1Q3 streaming operator with interleaved velocity organization
    
    Distribution function order:
    [f0_x1, f1_x1, f2_x1, f0_x2, f1_x2, f2_x2, ..., f0_xn, f1_xn, f2_xn]
    
    Args:
        nx: number of spatial grid points
        hx: grid spacing
        
    Returns:
        Single streaming matrix S of size (3*nx, 3*nx)
    """
    
    # D1Q3 velocity vectors
    e = [
        [0],    # e0: rest particle
        [1],    # e1: moving right (+x)  
        [-1]    # e2: moving left (-x)
    ]
    
    n_total = 3 * nx  # Total degrees of freedom
    
    # Sparse matrix construction arrays
    I_idx = Int[]
    J_idx = Int[]
    V_vals = Float64[]
    
    # Helper function: get global index for velocity i at position j
    global_index(vel_idx, pos_idx) = (pos_idx - 1) * 3 + vel_idx
    
    # Construct streaming matrix
    for pos in 1:nx  # Loop over spatial positions
        for vel in 1:3  # Loop over velocities at each position
            
            row_idx = global_index(vel, pos)
            ex = e[vel][1]  # x-component of velocity
            
            if ex == 0
                # Rest particle: f0 doesn't stream, just stays
                push!(I_idx, row_idx)
                push!(J_idx, row_idx) 
                push!(V_vals, 0.0)
                
            else
                # Streaming particle: apply central difference
                # ∂f/∂x ≈ (f[pos+1] - f[pos-1])/(2*hx)
                
                if pos > 1  # Left neighbor exists
                    col_idx = global_index(vel, pos - 1)
                    push!(I_idx, row_idx)
                    push!(J_idx, col_idx)
                    push!(V_vals, -ex / (2 * hx))
                end
                
                if pos < nx  # Right neighbor exists
                    col_idx = global_index(vel, pos + 1)
                    push!(I_idx, row_idx)
                    push!(J_idx, col_idx)
                    push!(V_vals, ex / (2 * hx))
                end
            end
        end
    end
    
    S = sparse(I_idx, J_idx, V_vals, n_total, n_total)
    return S, e
end

# Example usage
nx = 5
hx = 1
S_interleaved, e_velocities = streaming_operator_D1Q3_interleaved(nx, hx)

println("D1Q3 Interleaved Streaming Operator")
println("=" ^ 50)
println("Grid points: $nx")
println("Grid spacing: $hx") 
println("Matrix size: $(size(S_interleaved))")
println("Organization: [f0_x1, f1_x1, f2_x1, f0_x2, f1_x2, f2_x2, ...]")
println()

println("Full streaming matrix:")
println(Matrix(S_interleaved))
println()

# Function to create the 1D LBE streaming matrix
function streaming_matrix_LBE(n::Int)
    # Define the base blocks for f2 (no movement), f3 (right movement), f1 (left movement)
    S2 = Matrix{Float64}(1I, n, n)  # Identity for f2 (no movement)

    # Right shift block with periodic boundary conditions
    S3 = circshift(S2, (0, -1))

    # Left shift block with periodic boundary conditions
    S1 = circshift(S2, (0, 1))
    #= To match the definition of phi(x), return S' =  
    [f1(x1), f2(x1), ..., fQ(x1)];
    [f1(x2), f2(x2), ..., fQ(x2)];
    .
    [f1(xn), f2(xn), ..., fQ(xn)];
    =# 
    # Combine the blocks into the full streaming matrix
    S = 
         [zeros(n, n)  zeros(n, n) S1;
         S2  zeros(n, n)  zeros(n, n);
         zeros(n, n)  S3  zeros(n, n)]
         
    return S'
end

