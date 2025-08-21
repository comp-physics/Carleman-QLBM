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

function streaming_operator_D2Q9_interleaved(nx, ny, hx, hy)
    """
    Create D2Q9 streaming operator with interleaved velocity organization
    
    Distribution function order:
    [f0_x1y1, f1_x1y1, f2_x1y1, ..., f8_x1y1, f0_x2y1, f1_x2y1, ..., f8_x2y1, 
     f0_x3y1, ..., f0_x1y2, f1_x1y2, ..., f8_xnyn]
    
    Args:
        nx, ny: number of spatial grid points in x, y directions
        hx, hy: grid spacing in x, y directions
        
    Returns:
        Single streaming matrix S of size (9*nx*ny, 9*nx*ny)
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
    
    n_velocities = 9
    n_spatial = nx * ny
    n_total = n_velocities * n_spatial
    
    # Sparse matrix construction arrays
    I_idx = Int[]
    J_idx = Int[]
    V_vals = Float64[]
    
    # Helper function: get global index for velocity vel at position (i,j)
    global_index(vel_idx, i, j) = ((j-1) * nx + (i-1)) * n_velocities + vel_idx
    
    # Construct streaming matrix
    for j in 1:ny      # Loop over y positions
        for i in 1:nx  # Loop over x positions
            for vel in 1:n_velocities  # Loop over velocities at each position
                
                row_idx = global_index(vel, i, j)
                ex, ey = e[vel]  # x and y components of velocity
                
                if ex == 0 && ey == 0
                    # Rest particle: f0 doesn't stream
                    push!(I_idx, row_idx)
                    push!(J_idx, row_idx) 
                    push!(V_vals, 0.0)
                    
                else
                    # Streaming particle: apply central difference
                    # e⋅∇f = ex * ∂f/∂x + ey * ∂f/∂y
                    
                    # x-direction streaming: ex * ∂f/∂x
                    if ex != 0
                        if i > 1  # Left neighbor exists
                            col_idx = global_index(vel, i-1, j)
                            push!(I_idx, row_idx)
                            push!(J_idx, col_idx)
                            push!(V_vals, -ex / (2 * hx))
                        end
                        
                        if i < nx  # Right neighbor exists
                            col_idx = global_index(vel, i+1, j)
                            push!(I_idx, row_idx)
                            push!(J_idx, col_idx)
                            push!(V_vals, ex / (2 * hx))
                        end
                    end
                    
                    # y-direction streaming: ey * ∂f/∂y
                    if ey != 0
                        if j > 1  # Bottom neighbor exists
                            col_idx = global_index(vel, i, j-1)
                            push!(I_idx, row_idx)
                            push!(J_idx, col_idx)
                            push!(V_vals, -ey / (2 * hy))
                        end
                        
                        if j < ny  # Top neighbor exists
                            col_idx = global_index(vel, i, j+1)
                            push!(I_idx, row_idx)
                            push!(J_idx, col_idx)
                            push!(V_vals, ey / (2 * hy))
                        end
                    end
                end
            end
        end
    end
    
    S = sparse(I_idx, J_idx, V_vals, n_total, n_total)
    return S, e
end

# Example usage
nx, ny = 4, 3
hx, hy = 0.1, 0.1
S_interleaved_2D, e_velocities_2D = streaming_operator_D2Q9_interleaved(nx, ny, hx, hy)

println("D2Q9 Interleaved Streaming Operator")
println("=" ^ 50)
println("Grid size: $nx × $ny")
println("Grid spacing: hx=$hx, hy=$hy") 
println("Matrix size: $(size(S_interleaved_2D))")
println("Non-zeros: $(nnz(S_interleaved_2D))")
println("Organization: [f0_x1y1, f1_x1y1, ..., f8_x1y1, f0_x2y1, ...]")
println()

function visualize_D2Q9_structure(nx, ny)
    """
    Visualize the indexing structure for D2Q9 interleaved format
    """
    
    println("D2Q9 Interleaved Indexing Structure")
    println("=" ^ 45)
    println("Grid: $nx × $ny")
    println()
    
    # Show first few global indices
    println("Global index mapping (first 3×3 positions):")
    println("Position (i,j) → Global indices for [f0, f1, f2, ..., f8]")
    println("-" ^ 60)
    
    for j in 1:min(3, ny)
        for i in 1:min(3, nx)
            indices = [(((j-1) * nx + (i-1)) * 9 + vel) for vel in 1:9]
            println("($i,$j) → $(indices)")
        end
    end
    
    if nx > 3 || ny > 3
        println("...")
    end
    
    # Show velocity vectors
    println("\nVelocity vectors:")
    e = [[0,0], [1,0], [0,1], [-1,0], [0,-1], [1,1], [-1,1], [-1,-1], [1,-1]]
    for (i, vel) in enumerate(e)
        println("e$(i-1) = $vel")
    end
end

visualize_D2Q9_structure(nx, ny)

function test_D2Q9_interleaved_streaming()
    """
    Test the D2Q9 interleaved streaming operator with known functions
    """
    
    nx, ny = 5, 4
    hx, hy = 1.0, 1.0
    S, e = streaming_operator_D2Q9_interleaved(nx, ny, hx, hy)
    
    println("Testing D2Q9 Interleaved Streaming")
    println("=" ^ 40)
    println("Grid: $nx × $ny")
    
    # Create test distribution functions
    n_total = 9 * nx * ny
    f = zeros(n_total)
    
    # Initialize with different test functions for each velocity
    for j in 1:ny
        for i in 1:nx
            x_pos = Float64(i)
            y_pos = Float64(j)
            
            # Different test functions for each velocity
            for vel in 1:9
                global_idx = ((j-1) * nx + (i-1)) * 9 + vel
                
                if vel == 1      # f0: constant
                    f[global_idx] = 1.0
                elseif vel == 2  # f1: linear in x
                    f[global_idx] = x_pos
                elseif vel == 3  # f2: linear in y  
                    f[global_idx] = y_pos
                elseif vel == 4  # f3: x²
                    f[global_idx] = x_pos^2
                elseif vel == 5  # f4: y²
                    f[global_idx] = y_pos^2
                elseif vel == 6  # f5: xy
                    f[global_idx] = x_pos * y_pos
                else            # f6,f7,f8: more complex functions
                    f[global_idx] = sin(π * x_pos / nx) * cos(π * y_pos / ny)
                end
            end
        end
    end
    
    # Apply streaming
    f_streamed = S * f
    
    # Display results for a few positions
    println("\nResults at selected positions:")
    println("Position | Velocity | Initial | After Streaming | Expected")
    println("---------|----------|---------|-----------------|----------")
    
    test_positions = [(2,2), (3,2), (2,3)]
    for (i, j) in test_positions
        for vel in [1, 2, 3]  # Show first 3 velocities
            global_idx = ((j-1) * nx + (i-1)) * 9 + vel
            initial_val = f[global_idx]
            streamed_val = f_streamed[global_idx]
            
            # Calculate expected result analytically
            ex, ey = e[vel]
            expected = "..."
            if vel == 1  # Rest particle
                expected = "0"
            elseif vel == 2  # e1⋅∇(x) = 1⋅1 = 1
                expected = "1"
            elseif vel == 3  # e2⋅∇(y) = 1⋅1 = 1  
                expected = "1"
            end
            
            println("($i,$j)    | f$(vel-1)      | $(round(initial_val, digits=2))    | $(round(streamed_val, digits=3))        | $expected")
        end
        println("---------|----------|---------|-----------------|----------")
    end
    
    return f, f_streamed
end

f_initial_2D, f_result_2D = test_D2Q9_interleaved_streaming()

function extract_velocities_at_position_2D(f_interleaved, i, j, nx, ny)
    """
    Extract all 9 velocities at position (i,j)
    """
    start_idx = ((j-1) * nx + (i-1)) * 9 + 1
    return f_interleaved[start_idx:start_idx+8]
end

function extract_velocity_field_2D(f_interleaved, vel, nx, ny)
    """
    Extract specific velocity component at all positions as 2D array
    """
    result = zeros(nx, ny)
    for j in 1:ny
        for i in 1:nx
            global_idx = ((j-1) * nx + (i-1)) * 9 + vel
            result[i, j] = f_interleaved[global_idx]
        end
    end
    return result
end

function set_velocity_at_position_2D!(f_interleaved, i, j, vel, value, nx)
    """
    Set specific velocity component at position (i,j)
    """
    global_idx = ((j-1) * nx + (i-1)) * 9 + vel
    f_interleaved[global_idx] = value
end

function reshape_to_grid_format(f_interleaved, nx, ny)
    """
    Reshape interleaved format to [nx, ny, 9] array for visualization
    """
    result = zeros(nx, ny, 9)
    for j in 1:ny
        for i in 1:nx
            for vel in 1:9
                global_idx = ((j-1) * nx + (i-1)) * 9 + vel
                result[i, j, vel] = f_interleaved[global_idx]
            end
        end
    end
    return result
end

# Demonstration of utility functions
println("\nUtility Functions Demo")
println("=" ^ 25)

nx_demo, ny_demo = 3, 3
f_demo = collect(1.0:(9*nx_demo*ny_demo))  # Sequential numbering

println("Demo vector length: $(length(f_demo))")
println("Grid size: $nx_demo × $ny_demo")
println()

# Extract velocities at position (2,2)
velocities_22 = extract_velocities_at_position_2D(f_demo, 2, 2, nx_demo, ny_demo)
println("Velocities at position (2,2): $velocities_22")

# Extract velocity field for f1 (vel=2)
vel_field = extract_velocity_field_2D(f_demo, 2, nx_demo, ny_demo)
println("Velocity field f1:")
println(vel_field')

# Reshape to grid format
grid_format = reshape_to_grid_format(f_demo, nx_demo, ny_demo)
println("Shape of grid format: $(size(grid_format))")

function performance_analysis_D2Q9()
    """
    Analyze performance for different grid sizes
    """
    
    grid_sizes = [(10,10), (20,20), (30,30), (50,50)]
    
    println("D2Q9 Performance Analysis")
    println("=" ^ 30)
    println("Grid Size | Assembly Time | Matrix Size | Non-zeros | Memory (MB)")
    println("----------|---------------|-------------|-----------|------------")
    
    for (nx, ny) in grid_sizes
        hx, hy = 1.0/nx, 1.0/ny
        
        # Time assembly
        t = @elapsed S, e = streaming_operator_D2Q9_interleaved(nx, ny, hx, hy)
        
        # Calculate memory usage (rough estimate)
        memory_mb = nnz(S) * 16 / (1024^2)  # 16 bytes per entry (8 for int indices, 8 for float value)
        
        println("$(nx)×$(ny)     | $(round(t*1000, digits=1)) ms        | $(size(S)[1])×$(size(S)[2]) | $(nnz(S))     | $(round(memory_mb, digits=2))")
        
        if nx >= 50  # Don't test larger sizes in demo
            break
        end
    end
end

performance_analysis_D2Q9()


function streaming_operator_D3Q27_interleaved(nx, ny, nz, hx, hy, hz)
    """
    Create D3Q27 streaming operator with interleaved velocity organization
    
    Distribution function order:
    [f0_x1y1z1, f1_x1y1z1, f2_x1y1z1, ..., f26_x1y1z1, 
     f0_x2y1z1, f1_x2y1z1, ..., f26_x2y1z1,
     f0_x3y1z1, ..., 
     f0_x1y2z1, f1_x1y2z1, ..., f26_x1y2z1,
     ..., 
     f0_xnynzn, f1_xnynzn, ..., f26_xnynzn]
    
    Args:
        nx, ny, nz: number of spatial grid points in x, y, z directions
        hx, hy, hz: grid spacing in x, y, z directions
        
    Returns:
        Single streaming matrix S of size (27*nx*ny*nz, 27*nx*ny*nz)
    """
    
    # D3Q27 velocity vectors: all combinations of {-1, 0, 1}³
    e = []
    for ez in [-1, 0, 1]
        for ey in [-1, 0, 1]
            for ex in [-1, 0, 1]
                push!(e, [ex, ey, ez])
            end
        end
    end
    
    n_velocities = 27
    n_spatial = nx * ny * nz
    n_total = n_velocities * n_spatial
    
    # Sparse matrix construction arrays
    I_idx = Int[]
    J_idx = Int[]
    V_vals = Float64[]
    
    # Helper function: get global index for velocity vel at position (i,j,k)
    global_index(vel_idx, i, j, k) = (((k-1) * ny + (j-1)) * nx + (i-1)) * n_velocities + vel_idx
    
    println("Building D3Q27 streaming operator...")
    println("Grid: $nx × $ny × $nz")
    println("Total DOF: $n_total")
    
    # Construct streaming matrix
    for k in 1:nz      # Loop over z positions
        for j in 1:ny  # Loop over y positions
            for i in 1:nx  # Loop over x positions
                for vel in 1:n_velocities  # Loop over velocities at each position
                    
                    row_idx = global_index(vel, i, j, k)
                    ex, ey, ez = e[vel]  # x, y, z components of velocity
                    
                    if ex == 0 && ey == 0 && ez == 0
                        # Rest particle: f0 doesn't stream
                        push!(I_idx, row_idx)
                        push!(J_idx, row_idx) 
                        push!(V_vals, 0.0)
                        
                    else
                        # Streaming particle: apply central difference
                        # e⋅∇f = ex * ∂f/∂x + ey * ∂f/∂y + ez * ∂f/∂z
                        
                        # x-direction streaming: ex * ∂f/∂x
                        if ex != 0
                            if i > 1  # Left neighbor exists
                                col_idx = global_index(vel, i-1, j, k)
                                push!(I_idx, row_idx)
                                push!(J_idx, col_idx)
                                push!(V_vals, -ex / (2 * hx))
                            end
                            
                            if i < nx  # Right neighbor exists
                                col_idx = global_index(vel, i+1, j, k)
                                push!(I_idx, row_idx)
                                push!(J_idx, col_idx)
                                push!(V_vals, ex / (2 * hx))
                            end
                        end
                        
                        # y-direction streaming: ey * ∂f/∂y
                        if ey != 0
                            if j > 1  # Bottom neighbor exists
                                col_idx = global_index(vel, i, j-1, k)
                                push!(I_idx, row_idx)
                                push!(J_idx, col_idx)
                                push!(V_vals, -ey / (2 * hy))
                            end
                            
                            if j < ny  # Top neighbor exists
                                col_idx = global_index(vel, i, j+1, k)
                                push!(I_idx, row_idx)
                                push!(J_idx, col_idx)
                                push!(V_vals, ey / (2 * hy))
                            end
                        end
                        
                        # z-direction streaming: ez * ∂f/∂z
                        if ez != 0
                            if k > 1  # Back neighbor exists
                                col_idx = global_index(vel, i, j, k-1)
                                push!(I_idx, row_idx)
                                push!(J_idx, col_idx)
                                push!(V_vals, -ez / (2 * hz))
                            end
                            
                            if k < nz  # Front neighbor exists
                                col_idx = global_index(vel, i, j, k+1)
                                push!(I_idx, row_idx)
                                push!(J_idx, col_idx)
                                push!(V_vals, ez / (2 * hz))
                            end
                        end
                    end
                end
            end
        end
        
        # Progress indicator for large grids
        if nz > 5 && k % max(1, nz÷10) == 0
            println("  Progress: $(round(100*k/nz, digits=1))%")
        end
    end
    
    println("Assembling sparse matrix...")
    S = sparse(I_idx, J_idx, V_vals, n_total, n_total)
    println("✓ Assembly complete!")
    
    return S, e
end

# Example usage (small grid for demonstration)
nx, ny, nz = 3, 3, 3
hx, hy, hz = 0.1, 0.1, 0.1
S_interleaved_3D, e_velocities_3D = streaming_operator_D3Q27_interleaved(nx, ny, nz, hx, hy, hz)

println("\nD3Q27 Interleaved Streaming Operator")
println("=" ^ 50)
println("Grid size: $nx × $ny × $nz")
println("Grid spacing: hx=$hx, hy=$hy, hz=$hz") 
println("Matrix size: $(size(S_interleaved_3D))")
println("Non-zeros: $(nnz(S_interleaved_3D))")
println("Sparsity: $(round(100*nnz(S_interleaved_3D)/prod(size(S_interleaved_3D)), digits=3))%")
println("Organization: [f0_x1y1z1, f1_x1y1z1, ..., f26_x1y1z1, f0_x2y1z1, ...]")


function visualize_D3Q27_structure(nx, ny, nz)
    """
    Visualize the indexing structure for D3Q27 interleaved format
    """
    
    println("D3Q27 Interleaved Indexing Structure")
    println("=" ^ 45)
    println("Grid: $nx × $ny × $nz")
    println()
    
    # Show D3Q27 velocity vectors (first 10 for brevity)
    println("D3Q27 velocity vectors (showing first 10):")
    e = []
    for ez in [-1, 0, 1]
        for ey in [-1, 0, 1]
            for ex in [-1, 0, 1]
                push!(e, [ex, ey, ez])
            end
        end
    end
    
    for i in 1:min(10, 27)
        println("e$(i-1) = $(e[i])")
    end
    if length(e) > 10
        println("... (17 more velocities)")
    end
    
    # Show indexing for first few positions
    println("\nGlobal index mapping (first few positions):")
    println("Position (i,j,k) → Range of global indices")
    println("-" ^ 50)
    
    for k in 1:min(2, nz)
        for j in 1:min(2, ny)
            for i in 1:min(2, nx)
                start_idx = (((k-1) * ny + (j-1)) * nx + (i-1)) * 27 + 1
                end_idx = start_idx + 26
                println("($i,$j,$k) → [$start_idx:$end_idx]")
            end
        end
    end
    
    if nx > 2 || ny > 2 || nz > 2
        println("...")
    end
    
    # Memory estimation
    total_dof = 27 * nx * ny * nz
    matrix_elements = total_dof^2
    sparse_elements = 27 * nx * ny * nz * 7  # Rough estimate: 7 connections per velocity on average
    memory_dense_gb = matrix_elements * 8 / (1024^3)
    memory_sparse_mb = sparse_elements * 16 / (1024^2)  # 16 bytes per sparse entry
    
    println("\nMemory requirements:")
    println("Total DOF: $(total_dof)")
    println("Dense matrix: $(round(memory_dense_gb, digits=2)) GB")
    println("Sparse matrix (estimated): $(round(memory_sparse_mb, digits=1)) MB")
end

visualize_D3Q27_structure(nx, ny, nz)

function extract_velocities_at_position_3D(f_interleaved, i, j, k, nx, ny, nz)
    """
    Extract all 27 velocities at position (i,j,k)
    """
    start_idx = (((k-1) * ny + (j-1)) * nx + (i-1)) * 27 + 1
    return f_interleaved[start_idx:start_idx+26]
end

function extract_velocity_field_3D(f_interleaved, vel, nx, ny, nz)
    """
    Extract specific velocity component at all positions as 3D array
    """
    result = zeros(nx, ny, nz)
    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                global_idx = (((k-1) * ny + (j-1)) * nx + (i-1)) * 27 + vel
                result[i, j, k] = f_interleaved[global_idx]
            end
        end
    end
    return result
end

function set_velocity_at_position_3D!(f_interleaved, i, j, k, vel, value, nx, ny)
    """
    Set specific velocity component at position (i,j,k)
    """
    global_idx = (((k-1) * ny + (j-1)) * nx + (i-1)) * 27 + vel
    f_interleaved[global_idx] = value
end

function reshape_to_grid_format_3D(f_interleaved, nx, ny, nz)
    """
    Reshape interleaved format to [nx, ny, nz, 27] array for visualization
    """
    result = zeros(nx, ny, nz, 27)
    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                for vel in 1:27
                    global_idx = (((k-1) * ny + (j-1)) * nx + (i-1)) * 27 + vel
                    result[i, j, k, vel] = f_interleaved[global_idx]
                end
            end
        end
    end
    return result
end

function get_spatial_neighbors_3D(i, j, k, nx, ny, nz)
    """
    Get valid spatial neighbors for position (i,j,k)
    """
    neighbors = []
    
    # 6 face neighbors
    if i > 1; push!(neighbors, (i-1, j, k, "left")); end
    if i < nx; push!(neighbors, (i+1, j, k, "right")); end
    if j > 1; push!(neighbors, (i, j-1, k, "bottom")); end
    if j < ny; push!(neighbors, (i, j+1, k, "top")); end
    if k > 1; push!(neighbors, (i, j, k-1, "back")); end
    if k < nz; push!(neighbors, (i, j, k+1, "front")); end
    
    return neighbors
end

# Demonstration of utility functions
println("\nUtility Functions Demo")
println("=" ^ 25)

nx_demo, ny_demo, nz_demo = 2, 2, 2
total_dof = 27 * nx_demo * ny_demo * nz_demo
f_demo = collect(1.0:total_dof)  # Sequential numbering

println("Demo vector length: $(length(f_demo))")
println("Grid size: $nx_demo × $ny_demo × $nz_demo")
println()

# Extract velocities at position (1,1,1)
velocities_111 = extract_velocities_at_position_3D(f_demo, 1, 1, 1, nx_demo, ny_demo, nz_demo)
println("First few velocities at position (1,1,1): $(velocities_111[1:5])...")

# Extract velocity field for f1 (vel=2) 
vel_field_3d = extract_velocity_field_3D(f_demo, 2, nx_demo, ny_demo, nz_demo)
println("Velocity field f1 shape: $(size(vel_field_3d))")
println("f1 at all positions: $(vec(vel_field_3d))")

# Get neighbors
neighbors = get_spatial_neighbors_3D(1, 1, 1, nx_demo, ny_demo, nz_demo)
println("Neighbors of (1,1,1): $neighbors")

function test_D3Q27_interleaved_streaming()
    """
    Test the D3Q27 interleaved streaming operator with simple functions
    """
    
    nx, ny, nz = 4, 3, 3
    hx, hy, hz = 1.0, 1.0, 1.0
    
    println("Testing D3Q27 Interleaved Streaming")
    println("=" ^ 40)
    println("Grid: $nx × $ny × $nz")
    
    # Build streaming operator
    S, e = streaming_operator_D3Q27_interleaved(nx, ny, nz, hx, hy, hz)
    
    # Create test distribution functions
    n_total = 27 * nx * ny * nz
    f = zeros(n_total)
    
    # Initialize with simple test functions
    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                x_pos = Float64(i)
                y_pos = Float64(j)
                z_pos = Float64(k)
                
                for vel in 1:27
                    global_idx = (((k-1) * ny + (j-1)) * nx + (i-1)) * 27 + vel
                    
                    if vel == 1      # f0: constant
                        f[global_idx] = 1.0
                    elseif vel == 2  # f1: linear in x (e=[1,0,0])
                        f[global_idx] = x_pos
                    elseif vel == 3  # f2: linear in y (e=[0,1,0])
                        f[global_idx] = y_pos
                    elseif vel == 4  # f3: linear in z (e=[0,0,1])  
                        f[global_idx] = z_pos
                    elseif vel == 5  # f4: x² 
                        f[global_idx] = x_pos^2
                    elseif vel == 6  # f5: y²
                        f[global_idx] = y_pos^2
                    elseif vel == 7  # f6: z²
                        f[global_idx] = z_pos^2
                    else             # Other velocities: small perturbation
                        f[global_idx] = 0.1 * sin(π * x_pos / nx) * cos(π * y_pos / ny) * sin(π * z_pos / nz)
                    end
                end
            end
        end
    end
    
    # Apply streaming
    println("Applying streaming operator...")
    f_streamed = S * f
    
    # Analyze results at center position
    center_i, center_j, center_k = nx÷2 + 1, ny÷2 + 1, nz÷2 + 1
    
    println("\nResults at center position ($center_i,$center_j,$center_k):")
    println("Velocity | e_vector    | Initial | Streamed | Expected")
    println("---------|-------------|---------|----------|----------")
    
    for vel in 1:min(10, 27)  # Show first 10 velocities
        global_idx = (((center_k-1) * ny + (center_j-1)) * nx + (center_i-1)) * 27 + vel
        initial_val = f[global_idx]
        streamed_val = f_streamed[global_idx]
        
        ex, ey, ez = e[vel]
        expected = "..."
        
        if vel == 1  # Rest particle
            expected = "0"
        elseif vel == 2  # e[1,0,0]⋅∇(x) = 1
            expected = "1"
        elseif vel == 3  # e[0,1,0]⋅∇(y) = 1
            expected = "1"  
        elseif vel == 4  # e[0,0,1]⋅∇(z) = 1
            expected = "1"
        elseif vel == 5  # e[-1,0,0]⋅∇(x²) = -2x
            expected = "$(round(-2*center_i, digits=1))"
        end
        
        println("f$(vel-1)      | $([ex,ey,ez]) | $(round(initial_val, digits=2))    | $(round(streamed_val, digits=3))   | $expected")
    end
    
    return f, f_streamed, S
end

# Run test (with smaller grid to keep output manageable)
f_initial_3D, f_result_3D, S_test = test_D3Q27_interleaved_streaming()

function performance_analysis_D3Q27()
    """
    Analyze performance and memory usage for different grid sizes
    """
    
    # Use small grids due to large memory requirements
    grid_sizes = [(5,5,5), (8,8,8), (10,10,10)]
    
    println("D3Q27 Performance Analysis")
    println("=" ^ 35)
    println("Grid Size | Assembly Time | DOF      | Non-zeros | Memory (MB)")
    println("----------|---------------|----------|-----------|------------")
    
    for (nx, ny, nz) in grid_sizes
        if nx * ny * nz * 27 > 50000  # Skip if too large for demo
            println("$(nx)×$(ny)×$(nz)   | Skipped (too large for demo)")
            continue
        end
        
        hx, hy, hz = 1.0/nx, 1.0/ny, 1.0/nz
        
        # Time assembly
        println("Building $nx×$ny×$nz grid...")
        t = @elapsed S, e = streaming_operator_D3Q27_interleaved(nx, ny, nz, hx, hy, hz)
        
        # Calculate memory usage
        dof = 27 * nx * ny * nz
        memory_mb = nnz(S) * 16 / (1024^2)  # 16 bytes per sparse entry
        
        println("$(nx)×$(ny)×$(nz)   | $(round(t, digits=2)) s        | $(dof)    | $(nnz(S))   | $(round(memory_mb, digits=1))")
        
        # Clean up for memory
        S = nothing
        GC.gc()
    end
    
    # Theoretical scaling analysis
    println("\nTheoretical Scaling Analysis:")
    println("Grid Size | DOF        | Dense Matrix (GB) | Est. Sparse (MB)")
    println("----------|------------|-------------------|------------------")
    
    theoretical_sizes = [(10,10,10), (20,20,20), (50,50,50), (100,100,100)]
    for (nx, ny, nz) in theoretical_sizes
        dof = 27 * nx * ny * nz
        dense_gb = (dof^2 * 8) / (1024^3)
        sparse_mb = (dof * 7 * 16) / (1024^2)  # Assume ~7 connections per DOF
        
        println("$(nx)×$(ny)×$(nz)  | $(dof)   | $(round(dense_gb, digits=1))              | $(round(sparse_mb, digits=1))")
    end
end

performance_analysis_D3Q27()
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

