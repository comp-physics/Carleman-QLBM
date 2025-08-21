function lbm_cons()
    #=
    w0, w1, and w2 are weight factors for resting particle f0 (middle), left moving particle f1, and right moving particle, respectively  
    =#
    a = 1.; b = 3.; c = 9/2; d = -3/2;
    w_odd = 1/9; w_even = 1/36;
    #
    if LX == 1 || LY == 1
        D = 1;
        Q = 3;
        #
        w0 = 2/3; w1 = 1/6; w2 = 1/6
        w = [w0, w1, w2]
    else
        D = 2;
        Q = 9;
        #
        w0 = 4/9; w1 = w_odd; w3 = w_odd; w5 = w_even; w7 = w_even;
        w2 = w_odd; w4 = w_odd; w6 = w_even; w8 = w_even;
        w = [w0, w1, w2, w3, w4, w5, w6, w7, w8];
    end
    #

    if LX == 1 || LY == 1
        e = zeros(Q, 2);
        e[1,:] = [0 0]
        e[2,:] = [1 0]
        e[3,:] = [-1 0]
    else
        e = zeros(Q, D);
        e[1,:] = [0 0];
        e[2,:] = [1 0];
        e[3,:] = [-1 0];
        e[4,:] = [0 1];
        e[5,:] = [0 -1];
        e[6,:] = [1 1];
        e[7,:] = [1 -1];
        e[8,:] = [-1 1];
        e[9,:] = [-1 -1];
    end
    #
    return a, b, c, d, D, Q, w, e
end

function lbm_param()
  cs2 = 1/3
  return cs2
end
