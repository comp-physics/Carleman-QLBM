QCFD_SRC = ENV["QCFD_SRC"]  
QCFD_HOME = ENV["QCFD_HOME"]  
include(QCFD_SRC * "CLBM/timeMarching.jl")
include(QCFD_HOME * "/visualization/plot_CLBM_LBM.jl")
include(QCFD_SRC * "LBM/f_initial.jl")


function CLBM_collision_test(Q, omega, f, C, truncation_order, dt, tau_value, e_value, n_time, l_plot)
#---arbitrary initial condition
    println("Single point CLBM collision test of the Forets and Pouly Carleman linearization")
    u0 = 0.1
    f_ini = f_ini_test(u0)

    #
#    fT, _ = timeMarching_collision(omega, f, f_ini, tau_value, e_value, dt, n_time, l_plot)
    VT_f, VT, _, fT = timeMarching_collision_CLBM(omega, f, tau_value, Q, C, truncation_order, e_value, dt, f_ini, n_time, l_plot)
    #
    if l_plot
        close("all")
        figure(figsize=(10, 6)) 
        subplots_adjust(left = 0.1, right = 0.99, top = 0.95, bottom = 0.1, wspace = 0.35)
        plot_CLBM_LBM(fT, VT_f, n_time, "r", "CLBM", "LBM")
    end
    #
    return fT, VT_f, VT 
end

function CLBM_collision_test_sparse(Q, omega, f, truncation_order, dt, tau_value, e_value, n_time, l_plot)
#---arbitrary initial condition for sparse version
    println("Single point CLBM collision test using SPARSE Carleman matrix implementation")
    u0 = 0.1
    f_ini = f_ini_test(u0)

    # Use the optimized sparse time marching function
    VT_f, VT, _, fT = timeMarching_collision_CLBM_sparse(omega, f, tau_value, Q, truncation_order, e_value, dt, f_ini, n_time, l_plot)
    
    if l_plot
        close("all")
        figure(figsize=(10, 6)) 
        subplots_adjust(left = 0.1, right = 0.99, top = 0.95, bottom = 0.1, wspace = 0.35)
        plot_CLBM_LBM(fT, VT_f, n_time, "r", "CLBM (Sparse)", "LBM")
    end
    
    return fT, VT_f, VT 
end
