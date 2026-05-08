l_sympy = true
QCFD_SRC = ENV["QCFD_SRC"]
QCFD_HOME = ENV["QCFD_HOME"]

using PyPlot
using LaTeXStrings
using Printf
using LinearAlgebra
using SparseArrays
using Statistics

include(QCFD_SRC * "CLBE/clbe_tg2d_run.jl")
include(QCFD_SRC * "LBE/direct_LBE.jl")

"""
    main(; k_values=[3, 4], nx=3, ny=3, amplitude=0.02, local_n_time=100, coeff_method=:numerical, integrator=:euler)

D2Q9 mirror of `plot_truncation_order_error_comparison.jl`: compare the
domain-averaged absolute and relative error of each velocity component
f_m (m = 1..9) for several Carleman truncation orders k, CLBM vs direct
n-point LBE (apples-to-apples, same F, same S, same dt).

rho_value is pinned to 1.0001 (matching the D1Q3 multigrid strategy) so the
polynomial decomposition is well-posed and the ⟨ρ⟩ ≈ 1 Taylor assumption
of Li et al. holds. The `:matrix_exponential` integrator uses a sparse
Krylov expv-style propagator for large lifted systems.
"""
function main(; k_values=[3, 4], nx=3, ny=3, amplitude=0.02, local_n_time=100, coeff_method=:numerical, integrator=:euler, direct_lbe_integrator=:euler)
    if length(k_values) != 2
        error("This simplified script is intended for exactly two truncation orders, e.g. k_values=[3, 4].")
    end

    # QCFD convention: ngrid = LX * LY * LZ.
    global LX = nx; global LY = ny; global LZ = 1
    global ngrid = LX * LY * LZ
    ngrid_local = ngrid
    global Q_local = 9
    global Q = 9; global D = 2; global use_sparse = true
    global force_factor = 0.0
    global rho0 = 1.0001
    global lTaylor = true
    # Explicit-Euler stability: override LBM-unit default dt=1.0 from
    # clbm_config.jl. CLBM and direct n-point LBE share this dt.
    global dt = 0.1

    avg_abs_err_by_k = Dict{Int, Matrix{Float64}}()
    avg_rel_err_by_k = Dict{Int, Matrix{Float64}}()

    for k in k_values
        global truncation_order = k
        setup = build_carleman_setup(rho_value=rho0, nspatial=ngrid_local, method=coeff_method)
        global w_value  = setup.numeric_weights
        global e_value  = setup.numeric_velocities
        global F1_ngrid = setup.carleman_F1
        global F2_ngrid = setup.carleman_F2
        global F3_ngrid = setup.carleman_F3

        phi_ini = tg2d_initial_condition(nx, ny, amplitude, rho0,
            w_value, e_value, 1.0, 3.0, 9.0/2.0, -3.0/2.0)
        S_lbm, _ = streaming_operator_D2Q9_interleaved_periodic(nx, ny, 1.0, 1.0)

        phiT_lbe = timeMarching_direct_LBE_ngrid(phi_ini, dt, local_n_time,
            F1_ngrid, F2_ngrid, F3_ngrid; S_lbm=S_lbm, integrator=direct_lbe_integrator)
        phiT_clbm, _ = timeMarching_state_CLBM_sparse(
            setup.symbolic_collision, setup.symbolic_state, 1.0, Q, k,
            dt, phi_ini, local_n_time; S_lbm=S_lbm, nspatial=ngrid_local, integrator=integrator,
        )

        avg_lbe  = zeros(Q, local_n_time)
        avg_clbm = zeros(Q, local_n_time)
        for nt = 1:local_n_time
            avg_lbe[:,  nt] = vec(mean(reshape(phiT_lbe[:,  nt],  Q, ngrid_local), dims=2))
            avg_clbm[:, nt] = vec(mean(reshape(phiT_clbm[:, nt], Q, ngrid_local), dims=2))
        end
        abs_err = abs.(avg_clbm .- avg_lbe)
        rel_err = abs_err ./ max.(abs.(avg_lbe), eps(Float64))

        avg_abs_err_by_k[k] = abs_err
        avg_rel_err_by_k[k] = rel_err
    end

    time_axis = 1:local_n_time
    color_cycle = ["b", "r", "g", "m", "c", "k"]
    linestyle_cycle = ["-", "--", "-.", ":"]
    marker_cycle = ["o", "s", "^", "d", "v", ">", "<", "p", "h"]
    colors = Dict(k => color_cycle[mod1(i, length(color_cycle))] for (i, k) in enumerate(k_values))
    linestyles = Dict(k => linestyle_cycle[mod1(i, length(linestyle_cycle))] for (i, k) in enumerate(k_values))
    markers = Dict(k => marker_cycle[mod1(i, length(marker_cycle))] for (i, k) in enumerate(k_values))

    close("all")
    figure(figsize=(22, 7))
    for m = 1:Q
        subplot(2, Q, m)
        for k in k_values
            semilogy(
                time_axis,
                max.(avg_abs_err_by_k[k][m, :], 1e-300),
                linestyle=linestyles[k], color=colors[k], linewidth=1.6,
                marker=markers[k], markersize=3,
                markerfacecolor=colors[k], markeredgecolor=colors[k],
                label="k = $k",
            )
        end
        xlabel("Time step")
        ylabel(latexstring("|\\langle f_{$m} \\rangle^{\\mathrm{CLBM}} - \\langle f_{$m} \\rangle^{\\mathrm{LBE}}|"))
        if m == 1; legend(loc="best"); end

        subplot(2, Q, Q + m)
        for k in k_values
            semilogy(
                time_axis,
                max.(avg_rel_err_by_k[k][m, :], 1e-300),
                linestyle=linestyles[k], color=colors[k], linewidth=1.6,
                marker=markers[k], markersize=3,
                markerfacecolor=colors[k], markeredgecolor=colors[k],
                label="k = $k",
            )
        end
        xlabel("Time step")
        ylabel(latexstring("|\\Delta\\langle f_{$m} \\rangle| / |\\langle f_{$m} \\rangle^{\\mathrm{LBE}}|"))
    end

    suptitle("D2Q9 truncation-order error comparison, nx = ny = $nx, A = $amplitude, ρ₀ = $rho0")
    tight_layout(rect=(0, 0, 1, 0.95))

    output_dir = get(ENV, "QCFD_QCLBM_FIG_DIR", joinpath(homedir(), "Documents", "git-tex", "QC", "QCFD-QCLBM", "figs"))
    mkpath(output_dir)
    output_file = joinpath(output_dir, "plot_truncation_order_error_comparison_D2Q9_nx$(nx).pdf")
    savefig(output_file, bbox_inches="tight")
    # Also save a copy under data/ for convenience.
    data_copy = joinpath("data", "plot_truncation_order_error_comparison_D2Q9_nx$(nx).pdf")
    mkpath("data")
    savefig(data_copy, bbox_inches="tight")

    println("\nCompared truncation orders: ", collect(k_values))
    for k in k_values
        @printf "truncation order k = %d\n" k
        @printf "  overall max domain-averaged absolute error = %.3e\n" maximum(avg_abs_err_by_k[k])
        @printf "  overall max domain-averaged relative error = %.3e\n" maximum(avg_rel_err_by_k[k])
    end
    println("Saved D2Q9 truncation-order error comparison plot to:")
    println("  ", output_file)
    println("  ", data_copy)
    return output_file
end
