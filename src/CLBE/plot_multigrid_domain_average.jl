l_sympy = true
QCFD_SRC = ENV["QCFD_SRC"]
QCFD_HOME = ENV["QCFD_HOME"]

using PyPlot
using HDF5
using LaTeXStrings

include(QCFD_HOME * "/visualization/plot_kit.jl")
include(QCFD_HOME * "/visualization/plot_CLBE_LBM.jl")

if l_sympy
    using SymPy
    using LinearAlgebra
    include(QCFD_SRC * "CLBE/coeffs_poly.jl")
else
    using Symbolics
end

include("clbe_multigrid_run.jl")

function multigrid_initial_condition(comparison_ngrid; initial_condition=:legacy, u_ini=0.1)
    return d1q3_multigrid_initial_condition(comparison_ngrid; initial_condition=initial_condition, u_ini=u_ini)
end

function compute_domain_average_comparison(; local_n_time=max(n_time, 40), comparison_ngrid=3, local_truncation_order=truncation_order, coeff_method=coeff_generation_method, initial_condition=:legacy, u_ini=0.1, dt_override=nothing, integrator=:euler)
    core_cfg, case_cfg = build_d1q3_multigrid_configs(
        comparison_ngrid=comparison_ngrid,
        local_use_sparse=true,
        local_n_time=local_n_time,
        local_truncation_order=local_truncation_order,
        coeff_method=coeff_method,
        initial_condition=initial_condition,
        u_ini=u_ini,
        dt_override=dt_override,
    )
    phiT_lbe, phiT_clbm, VT, avg_abs_err, avg_rel_err = run_d1q3_multigrid(case_cfg, core_cfg; l_plot=false, integrator=integrator)

    avg_lbe = domain_average_distribution_history(phiT_lbe, Q, comparison_ngrid)
    avg_clbm = domain_average_distribution_history(phiT_clbm, Q, comparison_ngrid)

    return (
        phi_ini=multigrid_initial_condition(comparison_ngrid; initial_condition=initial_condition, u_ini=u_ini),
        phiT_lbe=phiT_lbe,
        phiT_clbm=phiT_clbm,
        VT=VT,
        avg_lbe=avg_lbe,
        avg_clbm=avg_clbm,
        avg_abs_err=avg_abs_err,
        avg_rel_err=avg_rel_err,
        comparison_ngrid=comparison_ngrid,
        local_n_time=local_n_time,
        local_truncation_order=local_truncation_order,
        initial_condition=initial_condition,
        u_ini=u_ini,
    )
end

function plot_domain_average_comparison!(comparison_data; figure_size=(12, 10))
    avg_lbe = comparison_data.avg_lbe
    avg_clbm = comparison_data.avg_clbm
    avg_abs_err = comparison_data.avg_abs_err
    avg_rel_err = comparison_data.avg_rel_err
    comparison_ngrid = comparison_data.comparison_ngrid
    local_n_time = comparison_data.local_n_time
    local_truncation_order = comparison_data.local_truncation_order
    initial_condition = comparison_data.initial_condition

    time = 1:local_n_time
    figure(figsize=figure_size)
    for m = 1:Q
        subplot(3, Q, m)
        plot(time, avg_lbe[m, :], "ok", label="LBM")
        plot(time, avg_clbm[m, :], "-", color="r", linewidth=1.8, label="CLBM")
        xlabel("Time step")
        ylabel(latexstring("\\langle f_{$m} \\rangle"))
        if m == 1
            legend(loc="best")
        end

        subplot(3, Q, Q + m)
        semilogy(time, avg_abs_err[m, :], "-", color="b", linewidth=1.8)
        xlabel("Time step")
        ylabel(latexstring("|\\langle f_{$m} \\rangle^{\\mathrm{CLBM}} - \\langle f_{$m} \\rangle^{\\mathrm{LBM}}|"))

        subplot(3, Q, 2 * Q + m)
        semilogy(time, avg_rel_err[m, :], "-", color="g", linewidth=1.8)
        xlabel("Time step")
        ylabel(latexstring("|\\langle f_{$m} \\rangle^{\\mathrm{CLBM}} - \\langle f_{$m} \\rangle^{\\mathrm{LBM}}| / \\langle f_{$m} \\rangle^{\\mathrm{LBM}}"))
    end
    suptitle("Domain-averaged CLBM vs LBM, ngrid = $comparison_ngrid, k = $local_truncation_order, IC = $(initial_condition)")
    tight_layout(rect=(0, 0, 1, 0.96))

    return gcf(), avg_abs_err, avg_rel_err
end

"""
    main(; local_n_time, comparison_ngrid, local_truncation_order, output_basename, coeff_method, initial_condition, u_ini, dt_override, integrator)

Run domain-averaged CLBM vs LBM comparison for D1Q3 multigrid.

Arguments:
    local_n_time: Number of time steps
    comparison_ngrid: Number of grid points
    local_truncation_order: Carleman truncation order (overrides global truncation_order)
    output_basename: Output file basename (optional)
    coeff_method: Coefficient generation method (optional)
    initial_condition: D1Q3 initial-condition selector (`:legacy` or `:sinusoidal`)
    u_ini: Velocity amplitude used by the sinusoidal initializer
    dt_override: Optional explicit time step; if omitted, uses the multigrid stability convention `tau_value / 10`
    integrator: CLBE time integrator (`:euler` or `:matrix_exponential`; the exponential option uses a sparse Krylov expv-style propagator for large lifted systems)

Example usage:
    main(local_n_time=100, comparison_ngrid=6, local_truncation_order=4, initial_condition=:sinusoidal, u_ini=0.1, dt_override=0.05, integrator=:euler)
"""
function main(; local_n_time=max(n_time, 40), comparison_ngrid=3, local_truncation_order=truncation_order, output_basename=nothing, coeff_method=coeff_generation_method, initial_condition=:legacy, u_ini=0.1, dt_override=nothing, integrator=:euler)
    comparison_data = compute_domain_average_comparison(
        local_n_time=local_n_time,
        comparison_ngrid=comparison_ngrid,
        local_truncation_order=local_truncation_order,
        coeff_method=coeff_method,
        initial_condition=initial_condition,
        u_ini=u_ini,
        dt_override=dt_override,
        integrator=integrator,
    )
    _, avg_abs_err, avg_rel_err = plot_domain_average_comparison!(comparison_data)

    if output_basename === nothing
        output_basename = "plot_multigrid_domain_average_D1Q3_$(initial_condition)_ngrid$(comparison_ngrid)_k$(local_truncation_order)_nt$(local_n_time)"
    end

    output_dir = get(ENV, "QCFD_QCLBM_FIG_DIR", joinpath(homedir(), "Documents", "git-tex", "QC", "QCFD-QCLBM", "figs"))
    mkpath(output_dir)
    output_file = joinpath(output_dir, output_basename * ".pdf")
    savefig(output_file)

    println("n_time used for CLBE/LBM comparison = ", local_n_time)
    println("Initial condition = ", initial_condition)
    println("Sinusoidal amplitude u_ini = ", u_ini)
    println("Time step dt = ", dt_override === nothing ? d1q3_multigrid_stable_dt(default_clbe_core_config()) : dt_override)
    println("Max domain-averaged absolute difference = ", maximum(avg_abs_err))
    println("Max domain-averaged relative difference = ", maximum(avg_rel_err))
    for m = 1:Q
        println("Max |Δ⟨f_$m⟩| = ", maximum(avg_abs_err[m, :]))
        println("Max |Δ⟨f_$m⟩| / max(|⟨f_$m⟩|, eps) = ", maximum(avg_rel_err[m, :]))
    end
    println("Saved figure to: ", output_file)

    display(gcf())
    show()

    return output_file, avg_abs_err, avg_rel_err
end
