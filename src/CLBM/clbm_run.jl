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

# symbolic calculation does not work for Kronecker product because it cannot distinguish f1f2 and f2f1

tau_value, n_time = 1. , 100
LX = 1; LY = 1; LZ = 1;  force_factor = 0.; dt_force_over_lbm = 1.
ngrid = 1
#tau_value, n_time = .51 , 40
#tau_value, n_time = 1.5 , 200000 # stable for larger than t>=200000
#tau_value = 1.4 # unstable for larger than t>=200000
#dt = 1 #2023-10-30/Xiaolong: dt must be one
dt = tau_value ./ 10 #2023-10-30/Xiaolong: dt must be one in LBM. This is just for testing the LBM-collision without streaming
lscale_kvector_tobox = false

include(QCFD_SRC * "CLBM/collision_sym.jl")
include(QCFD_SRC * "CLBM/carleman_transferA.jl")

include(QCFD_SRC * "CLBM/carleman_transferA_ngrid.jl")
include(QCFD_SRC * "CLBM/LBM_const_subs.jl")
#include(QCFD_SRC * "LBM/julia/julia_lib/matrix_kit.jl")
include(QCFD_SRC * "LBM/lbm_cons.jl")
include(QCFD_SRC * "LBM/lbm_const_sym.jl")
include(QCFD_SRC * "CLBM/CLBM_collision_test.jl")
include(QCFD_SRC * "LBM/forcing.jl")

Q = 3
D = 1

w, e, w_value, e_value = lbm_const_sym()
c = 1
#rho = 1
rho0 = 1.0001 # Any arbitrary flow
lTaylor = true
lorder2 = false
l_ini_feq = false



f, omega, u, rho = collision(Q, D, w, e, rho0, lTaylor, lorder2)

poly_order = 3 

#truncation_order = 4
truncation_order = 3


l_plot = true

##
#C = carleman_C(Q, truncation_order, poly_order, f, omega, tau_value)
#
global F1_ngrid, F2_ngrid, F3_ngrid = get_coeff_LBM_Fi_ngrid(poly_order, Q, f, omega, tau_value, ngrid)
C, bt, F0 = carleman_C(Q, truncation_order, poly_order, f, omega, tau_value, force_factor, w_value, e_value)
V = carleman_V(f, truncation_order)
##
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
