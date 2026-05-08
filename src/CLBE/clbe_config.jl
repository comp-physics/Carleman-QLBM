# CLBE Configuration Parameters
# Centralized configuration to avoid duplication across files

# Global simulation parameters
global tau_value = 1.0
global n_time = 10
# dt = 1 is the LBM lattice-time unit (Li et al., Eq. (eq:clbm)). It is the
# value required to match the discrete LBM streaming+collision iteration
# one-for-one. Individual test/benchmark drivers may override this globally
# to `tau_value / 10` (or similar) when running collision-only tests, or when
# explicit-Euler stability on the lifted Carleman operator demands a smaller
# step — the imaginary streaming eigenvalues make explicit Euler unstable at
# dt = 1 on multi-step runs. See Xiaolong 2023-10-30 comment in the reference
# repo `QCFD/src/CLBM/carleman_Forets_Pouly.jl` for the convention.
global dt = 1.0

# Domain parameters
global LX = 1
global LY = 1
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
global Q = 3
global D = 1

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

Base.@kwdef struct CLBECoreConfig
    tau_value::Float64 = 1.0
    n_time::Int = 10
    dt::Float64 = 1.0
    use_sparse::Bool = true
    force_factor::Float64 = 0.0
    dt_force_over_lbm::Float64 = 1.0
    Q::Int = 3
    D::Int = 1
    poly_order::Int = 3
    truncation_order::Int = 3
    lscale_kvector_tobox::Bool = false
    coeff_generation_method::Symbol = :numerical
    rho0::Float64 = 1.0001
    lTaylor::Bool = true
    lorder2::Bool = false
    l_ini_feq::Bool = false
    u0::Float64 = 0.1
end

Base.@kwdef struct D1Q3MultigridConfig
    LX::Int = 1
    LY::Int = 1
    LZ::Int = 1
    n_time::Int = 10
    use_sparse::Bool = true
    truncation_order::Int = 3
    coeff_generation_method::Symbol = :numerical
    initial_condition::Symbol = :legacy
    u_ini::Float64 = 0.1
    dt::Float64 = 1.0
end

config_ngrid(cfg::D1Q3MultigridConfig) = cfg.LX * cfg.LY * cfg.LZ

function default_clbe_core_config()
    return CLBECoreConfig(
        tau_value=tau_value,
        n_time=n_time,
        dt=dt,
        use_sparse=use_sparse,
        force_factor=force_factor,
        dt_force_over_lbm=dt_force_over_lbm,
        Q=Q,
        D=D,
        poly_order=poly_order,
        truncation_order=truncation_order,
        lscale_kvector_tobox=lscale_kvector_tobox,
        coeff_generation_method=coeff_generation_method,
        rho0=rho0,
        lTaylor=lTaylor,
        lorder2=lorder2,
        l_ini_feq=l_ini_feq,
        u0=u0,
    )
end

function d1q3_multigrid_stable_dt(core_cfg::CLBECoreConfig)
    return core_cfg.tau_value / 10
end

function default_d1q3_multigrid_config(; comparison_ngrid=3,
    local_n_time=max(n_time, 40),
    use_sparse_val=use_sparse,
    local_truncation_order=truncation_order,
    coeff_method=coeff_generation_method,
    initial_condition=:legacy,
    u_ini=0.1,
    dt_value=d1q3_multigrid_stable_dt(default_clbe_core_config()))

    return D1Q3MultigridConfig(
        LX=comparison_ngrid,
        LY=1,
        LZ=1,
        n_time=local_n_time,
        use_sparse=use_sparse_val,
        truncation_order=local_truncation_order,
        coeff_generation_method=coeff_method,
        initial_condition=initial_condition,
        u_ini=u_ini,
        dt=dt_value,
    )
end

function apply_clbe_core_config!(cfg::CLBECoreConfig)
    global tau_value = cfg.tau_value
    global n_time = cfg.n_time
    global dt = cfg.dt
    global use_sparse = cfg.use_sparse
    global force_factor = cfg.force_factor
    global dt_force_over_lbm = cfg.dt_force_over_lbm
    global Q = cfg.Q
    global D = cfg.D
    global poly_order = cfg.poly_order
    global truncation_order = cfg.truncation_order
    global lscale_kvector_tobox = cfg.lscale_kvector_tobox
    global coeff_generation_method = cfg.coeff_generation_method
    global rho0 = cfg.rho0
    global lTaylor = cfg.lTaylor
    global lorder2 = cfg.lorder2
    global l_ini_feq = cfg.l_ini_feq
    global u0 = cfg.u0
    return cfg
end

function apply_d1q3_multigrid_config!(cfg::D1Q3MultigridConfig)
    global LX = cfg.LX
    global LY = cfg.LY
    global LZ = cfg.LZ
    global ngrid = config_ngrid(cfg)
    global use_sparse = cfg.use_sparse
    global truncation_order = cfg.truncation_order
    global coeff_generation_method = cfg.coeff_generation_method
    global dt = cfg.dt
    return cfg
end

function configure_d1q3_runtime!(core_cfg::CLBECoreConfig, case_cfg::D1Q3MultigridConfig)
    apply_clbe_core_config!(core_cfg)
    apply_d1q3_multigrid_config!(case_cfg)
    return nothing
end

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
