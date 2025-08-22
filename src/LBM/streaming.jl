# LBD(:,:,:,2) is the f_m after collision, which is calculated from LBD(:,:,:,1). The streaming term is calculated by \vec{e_m}*(\partial f_m/ \partial (\vec{x_m})). Due to the periodic boundary conditions, the streaming is just a shift of LBD(:,:,:,2) according to \vec{e_m} since \vec{x_m} = \vec{e_m}*dt with dt=1. 
# Apply no-slip boundary conditions
function find_index_bounce_velocity(e, iQ)
    Q = size(e)[1]
    ex = e[:, 1]
    ey = e[:, 2]
    ey_k = ey[iQ]
    ex_k = ex[iQ]
    k_bounce = NaN
    for kk in 2:Q
        if ey[kk] == -ey_k && ex[kk] == -ex_k
            k_bounce = kk
        end
    end
    return k_bounce
end
#
function apply_no_slip_boundaries(LBD)
    ex = e[:, 1]
    ey = e[:, 2]
    for x in 1:LX
        for k in 1:Q
            # Top wall (no-slip bounce-back)
            if ey[k] == 1
               k_bounce = find_index_bounce_velocity(e, k)
               LBD[x, LY, k, 1] = LBD[x, LY, k_bounce, 2]
            end
            # Bottom wall (no-slip bounce-back)
            if ey[k] == -1
                k_bounce = find_index_bounce_velocity(e, k)
                LBD[x, 1, k, 1] = LBD[x, 1, k_bounce, 2]
            end
        end
    end
    return LBD 
end


function streaming(LBD, LX, LY, l_noslipBC)

    if l_noslipBC
        # 
        ind_top = LY - 1
        ind_bot = 2

        upindex = collect(ind_bot:ind_top+1)
        downindex = collect(ind_bot-1:ind_top)

        #
        LBD = apply_no_slip_boundaries(LBD)
        #
    else

        ind_top = LY
        ind_bot = 1

        upindex = push!(collect(ind_bot + 1:ind_top), 1)
        downindex = pushfirst!(collect(ind_bot:ind_top - 1), LY)

#        upindex = push!(collect(2:LY), 1)
#        downindex = pushfirst!(collect(1:LY-1), LY)
    end
    leftindex = pushfirst!(collect(1:LX-1), LX)
    rightindex = push!(collect(2:LX), 1);

    LBD[:,:,1,1] = LBD[:,:,1,2];
    LBD[:,:,2,1] = LBD[leftindex,:,2,2];
    LBD[:,:,3,1] = LBD[rightindex,:,3,2];
    if D > 1
#        LBD[:,:,4,1] = LBD[:,downindex,4,2];
#        LBD[:,:,5,1] = LBD[:,upindex,5,2];
#        LBD[:,:,6,1] = LBD[leftindex,downindex,6,2];
#        LBD[:,:,7,1] = LBD[leftindex,upindex,7,2];
#        LBD[:,:,8,1] = LBD[rightindex,downindex,8,2];
#        LBD[:,:,9,1] = LBD[rightindex,upindex,9,2];
#
        LBD[:,ind_bot:LY,4,1] = LBD[:,downindex,4,2];
        LBD[:,1:ind_top,5,1] = LBD[:,upindex,5,2];
        LBD[:,ind_bot:LY,6,1] = LBD[leftindex,downindex,6,2];
        LBD[:,1:ind_top,7,1] = LBD[leftindex,upindex,7,2];
        LBD[:,ind_bot:LY,8,1] = LBD[rightindex,downindex,8,2];
        LBD[:,1:ind_top,9,1] = LBD[rightindex,upindex,9,2];
        #
        # on in the middle
#        LBD[:,ind_bot:ind_top,4,1] = LBD[:,downindex,4,2];
#        LBD[:,ind_bot:ind_top,5,1] = LBD[:,upindex,5,2];
#        LBD[:,ind_bot:ind_top,6,1] = LBD[leftindex,downindex,6,2];
#        LBD[:,ind_bot:ind_top,7,1] = LBD[leftindex,upindex,7,2];
#        LBD[:,ind_bot:ind_top,8,1] = LBD[rightindex,downindex,8,2];
#        LBD[:,ind_bot:ind_top,9,1] = LBD[rightindex,upindex,9,2];
#        ##
#        LBD[:,ind_bot:ind_top,4,1] = LBD[:,upindex,4,2];
#        LBD[:,ind_bot:ind_top,5,1] = LBD[:,downindex,5,2];
#        LBD[:,ind_bot:ind_top,6,1] = LBD[leftindex,upindex,6,2];
#        LBD[:,ind_bot:ind_top,7,1] = LBD[leftindex,downindex,7,2];
#        LBD[:,ind_bot:ind_top,8,1] = LBD[rightindex,upindex,8,2];
#        LBD[:,ind_bot:ind_top,9,1] = LBD[rightindex,downindex,9,2];
        #
    end
    #
    streaming_term = LBD[:,:,:,1] .- LBD[:,:,:,2]
    #
    return LBD, streaming_term 
end
