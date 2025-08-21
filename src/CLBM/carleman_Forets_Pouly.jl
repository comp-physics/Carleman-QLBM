l_sympy = true
QCFD_SRC = ENV["QCFD_SRC"]  
QCFD_HOME = ENV["QCFD_HOME"]  

using HDF5
using PyPlot
using LaTeXStrings

include(QCFD_HOME * "/visualization/plot_kit.jl")

if l_sympy
    using SymPy
    using LinearAlgebra
    include(QCFD_SRC * "CLBM/coeffs_poly.jl")
else
    using Symbolics
end

function pcolor_matrix(L, nm_title)
    vmax = maximum(abs, L)
    vmin = - vmax
    figure(figsize=(8, 6.5))
    pcolormesh(L, cmap = "RdBu_r", shading = "auto", vmax = vmax, vmin = vmin)
    title(nm_title)
    colorbar()
    tight_layout()
end

function plot_eigval(C)
    figure(figsize=(7, 6))
    eigvalues = eigvals(C)
    println(maximum(maximum.(real(eigvalues))))
    plot(real(eigvalues), imag(eigvalues), ".k")
    ticklabel_format(style="sci", scilimits = (0,0), axis="x")
    xlabel(L"Re(\lambda)")
    ylabel(L"Im(\lambda)")
    tight_layout()
end

TEXPATH = ENV["TEXPATH"]
fig_dir = TEXPATH * "QC/CLBM_forced/fig/"

# symbolic calculation does not work for Kronecker product because it cannot distinguish f1f2 and f2f1

lscale_kvector_tobox = false
LX = 3; LY = 1; LZ = 1;  force_factor = 0.; dt_force_over_lbm = 1.

ngrid = LX * LY * LZ
#delta_alpha_ini = zeros(ngrid)   


include(QCFD_HOME * "/julia_lib/matrix_kit.jl")
include(QCFD_SRC * "LBM/lbm_cons.jl")
include(QCFD_SRC * "LBM/lbm_const_sym.jl")
include(QCFD_SRC * "LBM/forcing.jl")
include(QCFD_SRC * "CLBM/streaming_Carleman.jl")
include(QCFD_SRC * "CLBM/collision_sym.jl")
include(QCFD_SRC * "CLBM/carleman_transferA.jl")
include(QCFD_SRC * "CLBM/LBM_const_subs.jl")
include(QCFD_SRC * "CLBM/CLBM_collision_test.jl")
include(QCFD_SRC * "CLBM/carleman_transferA_ngrid.jl")

Q = 3
D = 1

w, e, w_value, e_value = lbm_const_sym()
c = 1
#rho = 1
rho0 = 1.0001 # Any arbitrary flow
lTaylor = true
lorder2 = false
l_ini_feq = false


tau_value, n_time = 1. , 100
#tau_value, n_time = .51 , 40
#tau_value, n_time = 1.5 , 200000 # stable for larger than t>=200000
#tau_value = 1.4 # unstable for larger than t>=200000
#dt = 1 #2023-10-30/Xiaolong: dt must be one
dt = tau_value ./ 10 #2023-10-30/Xiaolong: dt must be one in LBM. This is just for testing the LBM-collision without streaming

f, omega, u, rho = collision(Q, D, w, e, rho0, lTaylor, lorder2)

poly_order = 3 

truncation_order = 3

l_timeMarch = true
l_plot = true

F1_ngrid, F2_ngrid, F3_ngrid = get_coeff_LBM_Fi_ngrid(poly_order, Q, f, omega, tau_value, ngrid)

phi = get_phi(f, ngrid)
V_phi = carleman_V(phi, truncation_order) 

C, bt, F0 = carleman_C(Q, truncation_order, poly_order, f, omega, tau_value, force_factor, w_value, e_value)
V = carleman_V(f, truncation_order)
#

if ngrid > 1
    #S_LBM = streaming_matrix(ngrid)
    S_LBM, _ = streaming_operator_D1Q3_interleaved(ngrid, 1)
    S_Fj = get_S_Fj(S_LBM, ngrid)
    S = carleman_S(Q, truncation_order, poly_order, ngrid, S_Fj)
    C_full = C .- S 
else    
    C_full = C 
end

# eigval_C_full = real(eigvals(C_full))
# println("maximum(eigval_C_full)", maximum(eigval_C_full))
# t_arbitrary = collect(0:1.e5:1.e6)
# norm_exp_C = cal_matrix_exp(t_arbitrary, C_full)
#
if l_timeMarch
#---CLBM vs LBM---
    fT, VT_f, VT = CLBM_collision_test(Q, omega, f, C, truncation_order, dt, tau_value, e_value, n_time, l_plot)
            #
    title("CLBM-D1Q3, " * L"\tau=" *string(tau_value)  * L", u_0 = 0.1")

    lsavef = false

    if lsavef
        home = ENV["HOME"]
        dir_fig = home * "/Documents/git-tex/QC/CLBM_forced/fig/"
        fn_fig =  "CLBM_collision_D1Q3_tau" * string(tau_value) 
        savefig(dir_fig * fn_fig * ".png", dpi=300)
    end		
end
