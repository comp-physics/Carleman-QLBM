# CLBE Configuration Parameters
# Centralized configuration to avoid duplication across files

# Global simulation parameters
global tau_value = 1.0
global n_time = 600
# dt = 1 is the LBM lattice-time unit (Li et al., Eq. (eq:clbm)). It is the
# value required to match the discrete LBM streaming+collision iteration
# one-for-one. Individual test/benchmark drivers may override this globally
# to `tau_value / 10` (or similar) when running collision-only tests, or when
# explicit-Euler stability on the lifted Carleman operator demands a smaller
# step — the imaginary streaming eigenvalues make explicit Euler unstable at
# dt = 1 on multi-step runs. See Xiaolong 2023-10-30 comment in the reference
# repo `QCFD/src/CLBM/carleman_Forets_Pouly.jl` for the convention.
global dt = 0.1

# Domain parameters
global LX = 3
global LY = 3
global LZ = 1

# Grid parameters: derived from the physical lattice dimensions so that
# ngrid = LX * LY * LZ holds everywhere in the code (QCFD convention).
# Drivers that want a nonuniform ngrid should set LX/LY/LZ first, then
# re-derive ngrid = LX*LY*LZ.
global ngrid = LX * LY * LZ

# true: Use sparse matrices (recommended for ngrid > 1)
# false: Use dense matrices (only feasible for small problems)
global use_sparse = true 

# Helper function to get recommended sparse setting
function get_recommended_sparse_setting(ngrid_val)
    """
    Returns recommended sparse setting based on problem size.
    
    Args:
        ngrid_val: Grid size parameter
        
    Returns:
        bool: true if sparse is recommended, false if dense is optimal
    """
    if ngrid_val <= 1
        return false  # Dense is fine and equivalent for small problems
    else
        return true   # Sparse essential for larger problems
    end
end

# Helper function to validate sparse setting
function validate_sparse_setting(use_sparse_val, ngrid_val)
    """
    Validates sparse setting and provides warnings if needed.
    
    Args:
        use_sparse_val: Current sparse setting
        ngrid_val: Grid size parameter
    """
    recommended = get_recommended_sparse_setting(ngrid_val)
    
    if !use_sparse_val && ngrid_val > 1
        matrix_size = carleman_C_dim(3, 3, ngrid_val)  # Estimate with typical values
        memory_mb = matrix_size^2 * 8 / 1024^2  # 8 bytes per Float64
        
        if memory_mb > 100  # More than 100 MB
            println("⚠️  WARNING: Dense matrices will require ~$(round(memory_mb, digits=1)) MB")
            println("   Consider setting use_sparse = true for ngrid = $ngrid_val")
        end
    elseif use_sparse_val && ngrid_val == 1
        println("ℹ️  INFO: Sparse matrices provide minimal benefit for ngrid = 1")
        println("   Dense matrices are equivalent and may be easier to debug")
    end
end

# Forcing parameters
global force_factor = 0.0
global dt_force_over_lbm = 1.0

# LBM parameters
global Q = 9
global D = 2

# Carleman parameters
global poly_order = 3
global truncation_order = 3

# Other flags
global lscale_kvector_tobox = false
global coeff_generation_method = :numerical

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

println("CLBE configuration loaded:")
println("  tau_value = $tau_value")
println("  n_time = $n_time") 
println("  dt = $dt")
println("  Q = $Q, D = $D")
println("  truncation_order = $truncation_order")
println("  poly_order = $poly_order")
println("  force_factor = $force_factor")
println("  ngrid = $ngrid")
println("  coeff_generation_method = $coeff_generation_method")

if ngrid == 2
    println("⚠️  ngrid = 2 is a degenerate periodic centered-difference case.")
    println("   Use ngrid >= 3 for meaningful multigrid collision+streaming validation.")
end
