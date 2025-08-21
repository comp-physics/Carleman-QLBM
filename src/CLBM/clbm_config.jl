# CLBM Configuration Parameters
# Centralized configuration to avoid duplication across files

# Global simulation parameters
global tau_value = 1.0
global n_time = 10
global dt = tau_value / 10

# Domain parameters  
global LX = 1
global LY = 1
global LZ = 1

# Grid parameters
global ngrid = 1

# Forcing parameters
global force_factor = 0.0
global dt_force_over_lbm = 1.0

# LBM parameters
global Q = 3
global D = 1

# Carleman parameters
global poly_order = 3
global truncation_order = 3

# Other flags
global lscale_kvector_tobox = false

# Physical parameters
global rho0 = 1.0001  # Any arbitrary flow
global lTaylor = true
global lorder2 = false
global l_ini_feq = false

# Initial condition parameter
global u0 = 0.1

# These will be set by lbm_const_sym() - initialize with defaults
global w_value = [2/3, 1/6, 1/6]
global e_value = [0.0, 1.0, -1.0]

println("CLBM configuration loaded:")
println("  tau_value = $tau_value")
println("  n_time = $n_time") 
println("  dt = $dt")
println("  Q = $Q, D = $D")
println("  truncation_order = $truncation_order")
println("  poly_order = $poly_order")
println("  force_factor = $force_factor")
