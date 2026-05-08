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
    include(QCFD_SRC * "CLBE/coeffs_poly.jl")
else
    using Symbolics
end

include("clbe_config.jl")

include(QCFD_SRC * "CLBE/collision_sym.jl")
include(QCFD_SRC * "CLBE/carleman_transferA.jl")
include(QCFD_SRC * "CLBE/carleman_transferA_ngrid.jl")
include(QCFD_SRC * "CLBE/LBM_const_subs.jl")
include(QCFD_SRC * "LBM/lbm_cons.jl")
include(QCFD_SRC * "LBM/lbm_const_sym.jl")
include(QCFD_SRC * "CLBE/CLBE_collision_test.jl")
include(QCFD_SRC * "LBM/forcing.jl")
include(QCFD_SRC * "LBE/direct_LBE.jl")

function multigrid_initial_condition(comparison_ngrid; initial_condition=:legacy, u_ini=0.1)
    return d1q3_multigrid_initial_condition(comparison_ngrid; initial_condition=initial_condition, u_ini=u_ini)
end

function build_d1q3_multigrid_configs(; comparison_ngrid=ngrid,
    local_use_sparse=use_sparse,
    local_n_time=n_time,
    local_truncation_order=truncation_order,
    coeff_method=coeff_generation_method,
    initial_condition=:legacy,
    u_ini=0.1,
    dt_override=nothing)

    core_cfg = default_clbe_core_config()
    resolved_use_sparse = local_use_sparse
    if comparison_ngrid > 1 && !local_use_sparse
        println("ℹ️  Multigrid collision+streaming validation uses the sparse CLBM path; recording use_sparse=true in the multigrid case config")
        resolved_use_sparse = true
    end

    resolved_dt = dt_override === nothing ? d1q3_multigrid_stable_dt(core_cfg) : dt_override
    case_cfg = default_d1q3_multigrid_config(
        comparison_ngrid=comparison_ngrid,
        local_n_time=local_n_time,
        use_sparse_val=resolved_use_sparse,
        local_truncation_order=local_truncation_order,
        coeff_method=coeff_method,
        initial_condition=initial_condition,
        u_ini=u_ini,
        dt_value=resolved_dt,
    )
    return core_cfg, case_cfg
end

function prepare_d1q3_runtime(core_cfg::CLBECoreConfig, case_cfg::D1Q3MultigridConfig)
    configure_d1q3_runtime!(core_cfg, case_cfg)

    w, e, w_val, e_val = lbm_const_sym()
    global w_value = w_val
    global e_value = e_val

    f, omega, u, rho = collision(core_cfg.Q, core_cfg.D, w, e, core_cfg.rho0, core_cfg.lTaylor, core_cfg.lorder2)
    local_ngrid = config_ngrid(case_cfg)
    global F1_ngrid, F2_ngrid, F3_ngrid = get_coeff_LBM_Fi_ngrid(
        core_cfg.poly_order,
        core_cfg.Q,
        f,
        omega,
        core_cfg.tau_value,
        local_ngrid;
        method=case_cfg.coeff_generation_method,
        w_value_input=w_value,
        e_value_input=e_value,
        rho_value_input=core_cfg.rho0,
        lTaylor_input=core_cfg.lTaylor,
        D_input=core_cfg.D,
    )

    return (
        w=w,
        e=e,
        f=f,
        omega=omega,
        w_value=w_value,
        e_value=e_value,
        ngrid=local_ngrid,
    )
end

function run_legacy_collision_driver(runtime, core_cfg::CLBECoreConfig, case_cfg::D1Q3MultigridConfig; l_plot=true)
    local_ngrid = runtime.ngrid
    local_use_sparse = case_cfg.use_sparse
    local_n_time = case_cfg.n_time
    local_truncation_order = case_cfg.truncation_order

    if local_ngrid == 1 && local_use_sparse
        println("ℹ️  Using the original dense legacy collision test for ngrid = 1 to preserve single-point behavior")
        local_use_sparse = false
    end

    validate_sparse_setting(local_use_sparse, local_ngrid)

    if local_use_sparse
        println("Using SPARSE Carleman matrix implementation (use_sparse=$local_use_sparse, ngrid=$local_ngrid)")
        fT, VT_f, VT = CLBM_collision_test_sparse(core_cfg.Q, runtime.omega, runtime.f, local_truncation_order, case_cfg.dt, core_cfg.tau_value, runtime.e_value, local_n_time, l_plot)
    else
        println("Using DENSE Carleman matrix implementation (use_sparse=$local_use_sparse, ngrid=$local_ngrid)")
        C, bt, F0 = carleman_C(core_cfg.Q, local_truncation_order, core_cfg.poly_order, runtime.f, runtime.omega, core_cfg.tau_value, core_cfg.force_factor, runtime.w_value, runtime.e_value)
        fT, VT_f, VT = CLBM_collision_test(core_cfg.Q, runtime.omega, runtime.f, C, local_truncation_order, case_cfg.dt, core_cfg.tau_value, runtime.e_value, local_n_time, l_plot)
    end

    title("CLBM-D1Q3, τ=" * string(core_cfg.tau_value) * ", u_0 = 0.1")
    display(gcf())
    show()

    return fT, VT_f, VT
end

function run_singlepoint_affine_driver(runtime, core_cfg::CLBECoreConfig, case_cfg::D1Q3MultigridConfig; l_plot=true, integrator=:matrix_exponential, direct_lbe_integrator=:euler)
    local_n_time = case_cfg.n_time
    local_truncation_order = case_cfg.truncation_order

    phi_ini = multigrid_initial_condition(1; initial_condition=case_cfg.initial_condition, u_ini=case_cfg.u_ini)
    S_lbm = zeros(core_cfg.Q, core_cfg.Q)

    phiT_lbe = timeMarching_direct_LBE_ngrid(
        phi_ini,
        case_cfg.dt,
        local_n_time,
        F1_ngrid,
        F2_ngrid,
        F3_ngrid;
        S_lbm=S_lbm,
        integrator=direct_lbe_integrator,
    )

    phiT_clbm, VT = timeMarching_state_CLBM_sparse(
        runtime.omega,
        runtime.f,
        core_cfg.tau_value,
        core_cfg.Q,
        local_truncation_order,
        case_cfg.dt,
        phi_ini,
        local_n_time;
        S_lbm=S_lbm,
        nspatial=1,
        integrator=integrator,
    )

    avg_lbe = domain_average_distribution_history(phiT_lbe, core_cfg.Q, 1)
    avg_clbm = domain_average_distribution_history(phiT_clbm, core_cfg.Q, 1)
    avg_abs_err = abs.(avg_clbm .- avg_lbe)
    avg_rel_err = avg_abs_err ./ max.(abs.(avg_lbe), eps(Float64))

    if l_plot
        close("all")
        time = 1:local_n_time
        figure(figsize=(12, 10))
        for m = 1:core_cfg.Q
            subplot(3, core_cfg.Q, m)
            plot(time, avg_lbe[m, :], "ok", label="direct LBE")
            plot(time, avg_clbm[m, :], "-", color="r", linewidth=1.8, label="CLBE")
            xlabel("Time step")
            ylabel(latexstring("\\langle f_{$m} \\rangle"))
            if m == 1
                legend(loc="best")
            end

            subplot(3, core_cfg.Q, core_cfg.Q + m)
            semilogy(time, avg_abs_err[m, :], "-", color="b", linewidth=1.8)
            xlabel("Time step")
            ylabel(latexstring("|\\langle f_{$m} \\rangle^{\\mathrm{CLBE}} - \\langle f_{$m} \\rangle^{\\mathrm{LBE}}|"))

            subplot(3, core_cfg.Q, 2 * core_cfg.Q + m)
            semilogy(time, avg_rel_err[m, :], "-", color="g", linewidth=1.8)
            xlabel("Time step")
            ylabel(latexstring("|\\langle f_{$m} \\rangle^{\\mathrm{CLBE}} - \\langle f_{$m} \\rangle^{\\mathrm{LBE}}| / \\max(|\\langle f_{$m} \\rangle^{\\mathrm{LBE}}|, \\varepsilon)"))
        end
        suptitle("single-point affine CLBE comparison, k = $local_truncation_order, IC = $(case_cfg.initial_condition), integrator = $(integrator)")
        tight_layout(rect=(0, 0, 1, 0.96))
        display(gcf())
        show()
    end

    println("Running single-point affine CLBE/direct-LBE comparison (ngrid = 1, integrator = $(integrator))")
    println("Initial condition = ", case_cfg.initial_condition)
    println("Sinusoidal amplitude u_ini = ", case_cfg.u_ini)
    println("Max domain-averaged absolute difference = ", maximum(avg_abs_err))
    println("Max domain-averaged relative difference = ", maximum(avg_rel_err))

    return phiT_lbe, phiT_clbm, VT, avg_abs_err, avg_rel_err
end

function d1q3_interleaved_site_to_exact_lbm(site_state)
    return Float64[site_state[2], site_state[3], site_state[1]]
end

function d1q3_exact_lbm_site_to_interleaved(site_state)
    return Float64[site_state[3], site_state[1], site_state[2]]
end

function d1q3_exact_lbm_history(phi_ini, local_ngrid, local_n_time; tau_value)
    q_local = 3
    weights = (2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0)
    current_state = zeros(Float64, local_ngrid, q_local)
    collided_state = similar(current_state)
    streamed_state = similar(current_state)
    phi_hist = zeros(Float64, length(phi_ini), local_n_time)
    phi_hist[:, 1] .= Float64.(phi_ini)

    for ix = 1:local_ngrid
        block = phi_ini[(q_local * (ix - 1) + 1):(q_local * ix)]
        current_state[ix, :] .= d1q3_interleaved_site_to_exact_lbm(block)
    end

    for nt = 2:local_n_time
        for ix = 1:local_ngrid
            f0, fp, fm = current_state[ix, 1], current_state[ix, 2], current_state[ix, 3]
            rho = f0 + fp + fm
            u = abs(rho) <= eps(Float64) ? 0.0 : (fp - fm) / rho
            u2 = u^2

            feq0 = rho * weights[1] * (1.0 - 1.5 * u2)
            feqp = rho * weights[2] * (1.0 + 3.0 * u + 3.0 * u2)
            feqm = rho * weights[3] * (1.0 - 3.0 * u + 3.0 * u2)

            collided_state[ix, 1] = f0 - (f0 - feq0) / tau_value
            collided_state[ix, 2] = fp - (fp - feqp) / tau_value
            collided_state[ix, 3] = fm - (fm - feqm) / tau_value
        end

        for ix = 1:local_ngrid
            left_ix = ix == 1 ? local_ngrid : ix - 1
            right_ix = ix == local_ngrid ? 1 : ix + 1
            streamed_state[ix, 1] = collided_state[ix, 1]
            streamed_state[ix, 2] = collided_state[left_ix, 2]
            streamed_state[ix, 3] = collided_state[right_ix, 3]
        end

        current_state .= streamed_state

        for ix = 1:local_ngrid
            phi_hist[(q_local * (ix - 1) + 1):(q_local * ix), nt] .= d1q3_exact_lbm_site_to_interleaved(view(current_state, ix, :))
        end
    end

    avg_lbm_exact = domain_average_distribution_history(phi_hist, q_local, local_ngrid)
    return phi_hist, avg_lbm_exact
end

function run_multigrid_driver(runtime, core_cfg::CLBECoreConfig, case_cfg::D1Q3MultigridConfig; l_plot=true, integrator=:euler, direct_lbe_integrator=:euler, show_exact_lbm=false)
    local_ngrid = runtime.ngrid
    local_n_time = case_cfg.n_time
    local_truncation_order = case_cfg.truncation_order

    if local_ngrid == 2
        println("⚠️  ngrid = 2 is a degenerate periodic centered-difference case.")
        println("   The left and right neighbors coincide, so the centered streaming derivative collapses to zero.")
        println("   Use ngrid >= 3 for meaningful multigrid collision+streaming validation.")
    end

    validate_sparse_setting(case_cfg.use_sparse, local_ngrid)

    phi_ini = multigrid_initial_condition(local_ngrid; initial_condition=case_cfg.initial_condition, u_ini=case_cfg.u_ini)

    S_lbm, _ = streaming_operator_D1Q3_interleaved(local_ngrid, 1)
    phiT_lbe = timeMarching_direct_LBE_ngrid(phi_ini, case_cfg.dt, local_n_time, F1_ngrid, F2_ngrid, F3_ngrid; S_lbm=S_lbm, integrator=direct_lbe_integrator)
    phiT_clbm, VT = timeMarching_state_CLBM_sparse(runtime.omega, runtime.f, core_cfg.tau_value, core_cfg.Q, local_truncation_order, case_cfg.dt, phi_ini, local_n_time; S_lbm=S_lbm, nspatial=local_ngrid, integrator=integrator)

    avg_lbe = domain_average_distribution_history(phiT_lbe, core_cfg.Q, local_ngrid)
    avg_clbm = domain_average_distribution_history(phiT_clbm, core_cfg.Q, local_ngrid)
    avg_abs_err = abs.(avg_clbm .- avg_lbe)
    avg_rel_err = avg_abs_err ./ max.(abs.(avg_lbe), eps(Float64))
    show_exact_lbm_overlay = show_exact_lbm
    avg_lbm_exact = nothing

    if show_exact_lbm_overlay
        if isapprox(case_cfg.dt, 1.0; atol=1e-12, rtol=0.0)
            _, avg_lbm_exact = d1q3_exact_lbm_history(phi_ini, local_ngrid, local_n_time; tau_value=core_cfg.tau_value)
        else
            println("ℹ️  Exact discrete LBM comparison is only available for dt = 1; skipping exact-LBM overlay because dt = $(case_cfg.dt)")
            show_exact_lbm_overlay = false
        end
    end

    if l_plot
        close("all")
        time = 1:local_n_time
        figure(figsize=(12, 10))
        for m = 1:core_cfg.Q
            subplot(3, core_cfg.Q, m)
            if show_exact_lbm_overlay
                plot(time, avg_lbm_exact[m, :], "-k", linewidth=1.8, label="LBM")
                plot(time, avg_lbe[m, :], "or", markersize=5, markerfacecolor="none", label="LBE")
                plot(time, avg_clbm[m, :], "sb", markersize=5, markerfacecolor="none", label="CLBE")
            else
                plot(time, avg_lbe[m, :], "or", markersize=5, markerfacecolor="none", label="LBE")
                plot(time, avg_clbm[m, :], "-", color="b", linewidth=1.8, label="CLBE")
            end
            xlabel("Time step")
            ylabel(latexstring("\\langle f_{$m} \\rangle"))
            if m == 1
                legend(loc="best")
            end

            subplot(3, core_cfg.Q, core_cfg.Q + m)
            semilogy(time, avg_abs_err[m, :], "-", color="b", linewidth=1.8)
            xlabel("Time step")
            ylabel(latexstring("|\\langle f_{$m} \\rangle^{\\mathrm{CLBE}} - \\langle f_{$m} \\rangle^{\\mathrm{LBE}}|"))

            subplot(3, core_cfg.Q, 2 * core_cfg.Q + m)
            semilogy(time, avg_rel_err[m, :], "-", color="g", linewidth=1.8)
            xlabel("Time step")
            ylabel(latexstring("|\\langle f_{$m} \\rangle^{\\mathrm{CLBE}} - \\langle f_{$m} \\rangle^{\\mathrm{LBE}}| / \\max(|\\langle f_{$m} \\rangle^{\\mathrm{LBE}}|, \\varepsilon)"))
        end
        plot_label = show_exact_lbm_overlay ? "LBM / LBE / CLBE" : "LBE / CLBE"
        suptitle("ngrid = $local_ngrid, k = $local_truncation_order, IC = $(case_cfg.initial_condition), $plot_label")
        tight_layout(rect=(0, 0, 1, 0.96))
        display(gcf())
        show()
    end

    println("n_time used for CLBE/LBM multigrid comparison = ", local_n_time)
    println("Initial condition = ", case_cfg.initial_condition)
    println("Sinusoidal amplitude u_ini = ", case_cfg.u_ini)
    println("Max domain-averaged absolute difference = ", maximum(avg_abs_err))
    println("Max domain-averaged relative difference = ", maximum(avg_rel_err))
    for m = 1:core_cfg.Q
        println("Max |Δ⟨f_$m⟩| = ", maximum(avg_abs_err[m, :]))
        println("Max |Δ⟨f_$m⟩| / max(|⟨f_$m⟩|, eps) = ", maximum(avg_rel_err[m, :]))
    end

    return phiT_lbe, phiT_clbm, VT, avg_abs_err, avg_rel_err
end

"""
    main(; comparison_ngrid, local_use_sparse, local_n_time, local_truncation_order, l_plot, coeff_method, initial_condition, u_ini, dt_override, integrator, direct_lbe_integrator, show_exact_lbm)

Run D1Q3 CLBE/LBM multigrid driver with explicit truncation order.

Arguments:
    comparison_ngrid: Number of grid points
    local_use_sparse: Use sparse Carleman matrix
    local_n_time: Number of time steps
    local_truncation_order: Carleman truncation order (overrides global truncation_order)
    l_plot: Plot results
    coeff_method: Coefficient generation method (optional)
    initial_condition: D1Q3 initial-condition selector (`:legacy` or `:sinusoidal`)
    u_ini: Velocity amplitude used by the sinusoidal initializer
    dt_override: Optional explicit time step for the D1Q3 run; if omitted, uses the multigrid stability convention `tau_value / 10`
    integrator: CLBE time integrator (`:euler` or `:matrix_exponential`; the exponential option uses a sparse Krylov expv-style propagator for large lifted systems)
    direct_lbe_integrator: direct n-point LBE integrator (`:euler` or `:exponential_euler`, aliases `:etd`/`:etd1`)
    show_exact_lbm: When `true` and `dt = 1`, overlay the exact discrete D1Q3 LBM history in the comparison plot

Example usage:
    main(comparison_ngrid=6, local_truncation_order=4, local_n_time=100, initial_condition=:sinusoidal, u_ini=0.1, dt_override=1.0, integrator=:matrix_exponential, direct_lbe_integrator=:exponential_euler, show_exact_lbm=true)
"""
function run_d1q3_multigrid(case_cfg::D1Q3MultigridConfig, core_cfg::CLBECoreConfig=default_clbe_core_config(); l_plot=true, integrator=:euler, direct_lbe_integrator=:euler, show_exact_lbm=false)
    runtime = prepare_d1q3_runtime(core_cfg, case_cfg)
    integrator_key = normalize_clbe_integrator(integrator)

    if runtime.ngrid == 1
        if integrator_key == :euler
            println("Running legacy single-point CLBM collision test (ngrid = 1)")
            return run_legacy_collision_driver(runtime, core_cfg, case_cfg; l_plot=l_plot)
        end

        if show_exact_lbm
            println("ℹ️  show_exact_lbm is currently only used by the ngrid >= 2 D1Q3 multigrid plotting path; ignoring it for ngrid = 1")
        end

        return run_singlepoint_affine_driver(runtime, core_cfg, case_cfg; l_plot=l_plot, integrator=integrator_key, direct_lbe_integrator=direct_lbe_integrator)
    end

    println("Running validated multigrid CLBM collision+streaming comparison (ngrid = $(runtime.ngrid), coeff_method = $(case_cfg.coeff_generation_method), initial_condition = $(case_cfg.initial_condition))")
    return run_multigrid_driver(runtime, core_cfg, case_cfg; l_plot=l_plot, integrator=integrator, direct_lbe_integrator=direct_lbe_integrator, show_exact_lbm=show_exact_lbm)
end

function main(; comparison_ngrid=ngrid, local_use_sparse=use_sparse, local_n_time=n_time, local_truncation_order=truncation_order, l_plot=true, coeff_method=coeff_generation_method, initial_condition=:legacy, u_ini=0.1, dt_override=nothing, integrator=:euler, direct_lbe_integrator=:euler, show_exact_lbm=false)
    effective_n_time = local_n_time
    core_cfg, case_cfg = build_d1q3_multigrid_configs(
        comparison_ngrid=comparison_ngrid,
        local_use_sparse=local_use_sparse,
        local_n_time=effective_n_time,
        local_truncation_order=local_truncation_order,
        coeff_method=coeff_method,
        initial_condition=initial_condition,
        u_ini=u_ini,
        dt_override=dt_override,
    )
    return run_d1q3_multigrid(case_cfg, core_cfg; l_plot=l_plot, integrator=integrator, direct_lbe_integrator=direct_lbe_integrator, show_exact_lbm=show_exact_lbm)
end
