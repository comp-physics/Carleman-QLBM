using LinearAlgebra

function array_nan(len)
  urms_t = zeros(len) 
  urms_t[:] .= NaN
  return urms_t
end
#
function range_step(s, e, step)
    nstep = Int(round(abs((e - s) / step))) + 1
    return collect(range(s, e, nstep))                                            
end 
# 
function Kron_kth(ff, k)
    fk = ff
    if k == 1
        fk = fk
    else    
        for i=1:k-1
           fk = kron(ff, fk)
        end
    end
    return fk
end
#
function eigval_max(dum)
    lambda_max = maximum(maximum.(real(eigvals(dum)))) 
    return lambda_max
end
#
function norm_axis(m, axis)
  if axis == 1
    return norm.(eachrow(m))  # Use eachrow for row-wise norm
  else
    return norm.(eachcol(m))
  end
end
#
function cal_matrix_exp(t_arbitrary, C)
    norm_exp_C = zeros(length(t_arbitrary))
    for i = 1:length(t_arbitrary)
        norm_exp_C[i] = norm(exp(C * t_arbitrary[i]), Inf)
    end
    return norm_exp_C
end
