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

# Load centralized configuration
include("clbm_config.jl")

include(QCFD_SRC * "CLBM/collision_sym.jl")
include(QCFD_SRC * "CLBM/carleman_transferA.jl")

include(QCFD_SRC * "CLBM/carleman_transferA_ngrid.jl")
include(QCFD_SRC * "CLBM/LBM_const_subs.jl")
#include(QCFD_SRC * "LBM/julia/julia_lib/matrix_kit.jl")
include(QCFD_SRC * "LBM/lbm_cons.jl")
include(QCFD_SRC * "LBM/lbm_const_sym.jl")
include(QCFD_SRC * "CLBM/CLBM_collision_test.jl")
include(QCFD_SRC * "LBM/forcing.jl")

# Set up LBM constants (updates global w_value, e_value)
w, e, w_val, e_val = lbm_const_sym()
global w_value = w_val
global e_value = e_val

c = 1

f, omega, u, rho = collision(Q, D, w, e, rho0, lTaylor, lorder2)

l_plot = true

# Choose dense or sparse implementation based on ngrid size
global F1_ngrid, F2_ngrid, F3_ngrid = get_coeff_LBM_Fi_ngrid(poly_order, Q, f, omega, tau_value, ngrid)

# Validate and provide guidance on sparse setting
validate_sparse_setting(use_sparse, ngrid)

# Choose implementation based on user configuration
if use_sparse
    println("Using SPARSE Carleman matrix implementation (use_sparse=$use_sparse, ngrid=$ngrid)")
    fT, VT_f, VT = CLBM_collision_test_sparse(Q, omega, f, truncation_order, dt, tau_value, e_value, n_time, l_plot)
else
    println("Using DENSE Carleman matrix implementation (use_sparse=$use_sparse, ngrid=$ngrid)")
    C, bt, F0 = carleman_C(Q, truncation_order, poly_order, f, omega, tau_value, force_factor, w_value, e_value)
    V = carleman_V(f, truncation_order)
    fT, VT_f, VT = CLBM_collision_test(Q, omega, f, C, truncation_order, dt, tau_value, e_value, n_time, l_plot)
end

title("CLBM-D1Q3, " * L"\tau=" *string(tau_value)  * L", u_0 = 0.1")

lsavef = false

if lsavef
    home = ENV["HOME"]
    dir_fig = home * "/Documents/git-tex/QC/CLBM_forced/fig/"
    fn_fig =  "CLBM_collision_D1Q3_tau" * string(tau_value) 
    savefig(dir_fig * fn_fig * ".png", dpi=300)
end		
