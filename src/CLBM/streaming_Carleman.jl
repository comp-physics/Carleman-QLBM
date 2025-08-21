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

