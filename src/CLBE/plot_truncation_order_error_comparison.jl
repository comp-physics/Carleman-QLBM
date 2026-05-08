l_sympy = true
QCFD_SRC = ENV["QCFD_SRC"]
QCFD_HOME = ENV["QCFD_HOME"]

using PyPlot
using LaTeXStrings

include("plot_multigrid_domain_average.jl")

"""
    main(; k_values=[3, 4], comparison_ngrid=3, local_n_time=100, coeff_method=coeff_generation_method, initial_condition=:legacy, u_ini=0.1, dt_override=nothing, integrator=:euler)

Compare D1Q3 truncation-order behavior by reusing the multigrid
CLBE/LBM comparison workflow for each requested truncation order.

Arguments:
    k_values: Array of Carleman truncation orders to compare (e.g., [3, 4])
    comparison_ngrid: Number of grid points
    local_n_time: Number of time steps
    coeff_method: Coefficient generation method (optional)
    initial_condition: D1Q3 initial-condition selector (`:legacy` or `:sinusoidal`)
    u_ini: Velocity amplitude used by the sinusoidal initializer
    dt_override: Optional explicit time step shared by the compared runs
    integrator: CLBE time integrator (`:euler` or `:matrix_exponential`; the exponential option uses a sparse Krylov expv-style propagator for large lifted systems)

Example usage:
    main(k_values=[3,4], comparison_ngrid=6, local_n_time=100, initial_condition=:sinusoidal, u_ini=0.1, dt_override=0.05, integrator=:euler)
"""
function main(; k_values=[3, 4], comparison_ngrid=3, local_n_time=100, coeff_method=coeff_generation_method, initial_condition=:legacy, u_ini=0.1, dt_override=nothing, integrator=:euler)
    if length(k_values) != 2
        error("This simplified script is intended for exactly two truncation orders, e.g. k_values=[3, 4].")
    end

    avg_abs_err_by_k = Dict{Int, Matrix{Float64}}()
    avg_rel_err_by_k = Dict{Int, Matrix{Float64}}()

    for k in k_values
        comparison_data = compute_domain_average_comparison(
            local_n_time=local_n_time,
            comparison_ngrid=comparison_ngrid,
            local_truncation_order=k,
            coeff_method=coeff_method,
            initial_condition=initial_condition,
            u_ini=u_ini,
            dt_override=dt_override,
            integrator=integrator,
        )
        avg_abs_err_by_k[k] = comparison_data.avg_abs_err
        avg_rel_err_by_k[k] = comparison_data.avg_rel_err
    end

    time = 1:local_n_time
    color_cycle = ["b", "r", "g", "m", "c", "k"]
    linestyle_cycle = ["-", "--", "-.", ":"]
    marker_cycle = ["o", "s", "^", "d", "v", ">", "<", "p", "h"]
    colors = Dict(k => color_cycle[mod1(i, length(color_cycle))] for (i, k) in enumerate(k_values))
    linestyles = Dict(k => linestyle_cycle[mod1(i, length(linestyle_cycle))] for (i, k) in enumerate(k_values))
    markers = Dict(k => marker_cycle[mod1(i, length(marker_cycle))] for (i, k) in enumerate(k_values))

    close("all")
    figure(figsize=(12, 7))
    for m = 1:Q
        subplot(2, Q, m)
        for k in k_values
            semilogy(
                time,
                avg_abs_err_by_k[k][m, :],
                linestyle=linestyles[k],
                color=colors[k],
                linewidth=1.8,
                marker=markers[k],
                markersize=4,
                markerfacecolor=colors[k],
                markeredgecolor=colors[k],
                label="k = $k",
            )
        end
        xlabel("Time step")
        ylabel(latexstring("|\\langle f_{$m} \\rangle^{\\mathrm{CLBM}} - \\langle f_{$m} \\rangle^{\\mathrm{LBM}}|"))
        if m == 1
            legend(loc="best")
        end

        subplot(2, Q, Q + m)
        for k in k_values
            semilogy(
                time,
                avg_rel_err_by_k[k][m, :],
                linestyle=linestyles[k],
                color=colors[k],
                linewidth=1.8,
                marker=markers[k],
                markersize=4,
                markerfacecolor=colors[k],
                markeredgecolor=colors[k],
                label="k = $k",
            )
        end
        xlabel("Time step")
        ylabel(latexstring("|\\langle f_{$m} \\rangle^{\\mathrm{CLBM}} - \\langle f_{$m} \\rangle^{\\mathrm{LBM}}| / \\langle f_{$m} \\rangle^{\\mathrm{LBM}}"))
    end

    suptitle("Truncation-order error comparison, ngrid = $comparison_ngrid, IC = $(initial_condition)")
    tight_layout(rect=(0, 0, 1, 0.95))

    output_dir = get(ENV, "QCFD_QCLBM_FIG_DIR", joinpath(homedir(), "Documents", "git-tex", "QC", "QCFD-QCLBM", "figs"))
    mkpath(output_dir)
    output_file = joinpath(output_dir, "plot_truncation_order_error_comparison_D1Q3_$(initial_condition).pdf")
    savefig(output_file)

    println("Compared truncation orders: ", collect(k_values))
    println("Initial condition = ", initial_condition)
    println("Sinusoidal amplitude u_ini = ", u_ini)
    println("Time step dt = ", dt_override === nothing ? d1q3_multigrid_stable_dt(default_clbe_core_config()) : dt_override)
    for k in k_values
        println("truncation order k = ", k)
        println("  overall max domain-averaged absolute error = ", maximum(avg_abs_err_by_k[k]))
        println("  overall max domain-averaged relative error = ", maximum(avg_rel_err_by_k[k]))
    end
    println("Truncation order error comparison plot saved to: ", output_file)
    display(gcf())
    show()
    return output_file
end

