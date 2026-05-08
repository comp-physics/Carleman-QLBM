using Printf
using LinearAlgebra
using SparseArrays
using Statistics
using PyPlot
using LaTeXStrings

ENV["QCFD_HOME"] = pwd()
ENV["QCFD_SRC"]  = pwd() * "/src/"

include("src/CLBE/clbe_tg2d_run.jl")

# D2Q9 CLBE vs direct n-point LBE (apples-to-apples):
#   * domain-averaged density and velocity magnitude over time (both)
#   * absolute error of domain-averaged density and velocity magnitude
#   * per-site profile absolute error
#   * density and velocity error norms

function make_plot(nx, ny, k, output_path; local_n_time=n_time)
    result = run_tg2d_clbe_comparison(
        nx=nx,
        ny=ny,
        amplitude=0.02,
        rho_value=1.0001,
        local_n_time=local_n_time,
        boundary_setup=false,
        coeff_method=:numerical,
        local_truncation_order=k,
        reference_model=:direct_lbe,
    )
    save_dir = joinpath("data", @sprintf("d2q9_apples_plot_nx%d_ny%d_k%d_nt%d", nx, ny, k, local_n_time))
    saved_paths = save_tg2d_comparison_hdf5(result; output_dir=save_dir, snapshot_every=max(1, local_n_time ÷ 10))

    diag = result.diagnostics

    t = dt .* (0:local_n_time-1)
    safelog(x) = max.(x, 1e-10)

    close("all")
    fig = figure(figsize=(10, 8))

    ax = subplot(2, 2, 1)
    plot(t, diag.rho_ref_mean,  "-k", linewidth=2.0, label="direct n-point LBE")
    plot(t, diag.rho_clbm_mean, "--r", linewidth=1.6, label="CLBE")
    y_mid = mean(diag.rho_ref_mean); y_range = max(maximum(abs.(diag.rho_ref_mean .- y_mid)), 1e-8) * 5 + 1e-6
    ax.set_ylim(y_mid - y_range, y_mid + y_range)
    xlabel(L"t"); ylabel(L"\langle \rho \rangle")
    legend(loc="best"); grid(true)

    subplot(2, 2, 2)
    semilogy(t, safelog(diag.u_ref_rms),  "-k", linewidth=2.0, label="direct n-point LBE")
    semilogy(t, safelog(diag.u_clbm_rms), "--r", linewidth=1.6, label="CLBE")
    xlabel(L"t"); ylabel(L"\sqrt{\langle |u|^2 \rangle}")
    legend(loc="best"); grid(true)

    subplot(2, 2, 3)
    semilogy(t, safelog(diag.rho_mean_abs_err), "-b", linewidth=1.8)
    xlabel(L"t"); ylabel(L"|\langle\rho\rangle_{\mathrm{CLBE}} - \langle\rho\rangle_{\mathrm{LBE}}|")
    grid(true)

    subplot(2, 2, 4)
    semilogy(t, safelog(diag.u_rms_abs_err), "-b", linewidth=1.8)
    xlabel(L"t"); ylabel(L"|\langle|u|\rangle_{\mathrm{CLBE}} - \langle|u|\rangle_{\mathrm{LBE}}|")
    grid(true)

    suptitle("D2Q9 CLBE vs direct n-point LBE \nnx=$(nx), k=$(k), A=0.02, ρ₀=1.0001, nt=$(local_n_time), dt=$(dt)")
    tight_layout(rect=(0, 0, 1, 0.96))
    savefig(output_path, bbox_inches="tight")
    println("Saved plot: $output_path")
    println("Saved HDF5 snapshots/time series to: $(saved_paths.output_dir)")
    println("Saved time series to: $(saved_paths.time_series_path)")

    # Print key diagnostics so we can verify the plot didn't silently misreport
    @printf "\nDiagnostic values at the first and last step:\n"
    @printf "  ⟨ρ⟩_LBE[1]   = %.8e   ⟨ρ⟩_CLBE[1]   = %.8e\n" diag.rho_ref_mean[1]   diag.rho_clbm_mean[1]
    @printf "  ⟨ρ⟩_LBE[end] = %.8e   ⟨ρ⟩_CLBE[end] = %.8e\n" diag.rho_ref_mean[end] diag.rho_clbm_mean[end]
    @printf "  ⟨|u|⟩_LBE[1]   = %.8e   ⟨|u|⟩_CLBE[1]   = %.8e\n" diag.u_ref_rms[1]   diag.u_clbm_rms[1]
    @printf "  ⟨|u|⟩_LBE[end] = %.8e   ⟨|u|⟩_CLBE[end] = %.8e\n" diag.u_ref_rms[end] diag.u_clbm_rms[end]
    @printf "  max |Δ⟨ρ⟩|          = %.3e\n" maximum(diag.rho_mean_abs_err)
    @printf "  max |Δ⟨|u|⟩|        = %.3e\n" maximum(diag.u_rms_abs_err)
    @printf "  max density err norm = %.3e\n" maximum(diag.density_error_norm)
    @printf "  max velocity err norm = %.3e\n" maximum(diag.velocity_error_norm)
    @printf "  max |Δφ|∞           = %.3e\n" maximum(diag.profile_abs_max)
end

# By default, use `n_time` from src/CLBE/clbe_config_2D.jl. To override it for a
# specific plot, call for example `make_plot(3, 3, 3, output_path; local_n_time=100)`.
output_path = "data/d2q9_apples_to_apples_nx3_k3_nt$(n_time).pdf"
make_plot(3, 3, 3, output_path, local_n_time=200)
