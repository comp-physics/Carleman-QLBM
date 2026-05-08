using Printf
using Statistics
using LinearAlgebra
using PyPlot
using LaTeXStrings
using HDF5

ENV["QCFD_HOME"] = pwd()
ENV["QCFD_SRC"] = pwd() * "/src/"

include("src/CLBE/clbe_tg2d_run.jl")

function integrator_label(integrator)
    key = normalize_direct_lbe_integrator(integrator)
    if key == :euler
        return "Euler"
    elseif key == :exponential_euler
        return "Exponential Euler"
    end
    return string(key)
end

function integrator_key_string(integrator)
    return String(normalize_direct_lbe_integrator(integrator))
end

function kinetic_energy_history(phiT, nx, ny, e_value)
    n_time_local = size(phiT, 2)
    ke_mean = zeros(n_time_local)
    rho_mean = zeros(n_time_local)

    for nt in 1:n_time_local
        rho, ux, uy = macroscopic_fields_from_state(phiT[:, nt], nx, ny, e_value)
        rho_mean[nt] = mean(rho)
        ke_mean[nt] = mean(0.5 .* rho .* (ux .^ 2 .+ uy .^ 2))
    end

    return rho_mean, ke_mean
end

function run_lbe_integrator_histories(; nx=3,
    ny=3,
    amplitude=0.02,
    rho_value=1.0001,
    local_n_time=n_time,
    boundary_setup=false,
    coeff_method=:numerical,
    local_truncation_order=truncation_order,
    integrators=[:euler, :exponential_euler])

    if isempty(integrators)
        error("integrators must be non-empty")
    end

    runtime = prepare_d2q9_carleman_runtime(
        nx=nx,
        ny=ny,
        rho_value=rho_value,
        coeff_method=coeff_method,
        local_truncation_order=local_truncation_order,
        boundary_setup=boundary_setup,
    )

    phi_ini = tg2d_initial_condition(
        nx,
        ny,
        amplitude,
        rho_value,
        runtime.w_value,
        runtime.e_value,
        runtime.setup.a_val,
        runtime.setup.b_val,
        runtime.setup.c_val,
        runtime.setup.d_val,
    )

    histories = Dict{Symbol, Matrix{Float64}}()
    diagnostics = Dict{Symbol, NamedTuple}()

    for integrator in integrators
        key = normalize_direct_lbe_integrator(integrator)
        phiT = timeMarching_direct_LBE_ngrid(
            phi_ini,
            dt,
            local_n_time,
            runtime.setup.carleman_F1,
            runtime.setup.carleman_F2,
            runtime.setup.carleman_F3;
            S_lbm=runtime.S_lbm,
            integrator=key,
        )
        rho_mean, ke_mean = kinetic_energy_history(phiT, nx, ny, runtime.e_value)
        histories[key] = phiT
        diagnostics[key] = (
            rho_mean=rho_mean,
            kinetic_energy_mean=ke_mean,
        )
    end

    return (
        nx=nx,
        ny=ny,
        amplitude=amplitude,
        rho_value=rho_value,
        local_n_time=local_n_time,
        boundary_setup=boundary_setup,
        coeff_method=coeff_method,
        local_truncation_order=local_truncation_order,
        integrators=[normalize_direct_lbe_integrator(integ) for integ in integrators],
        phi_ini=phi_ini,
        histories=histories,
        diagnostics=diagnostics,
        e_value=runtime.e_value,
    )
end

function default_output_dir(; nx, ny, local_n_time, local_truncation_order)
    return joinpath("data", @sprintf("d2q9_lbe_integrator_compare_nx%d_ny%d_k%d_nt%d", nx, ny, local_truncation_order, local_n_time))
end

function save_lbe_integrator_time_series_hdf5(result; output_dir)
    mkpath(output_dir)
    ts_path = joinpath(output_dir, "time_series.h5")
    t = dt .* collect(0:result.local_n_time - 1)

    h5open(ts_path, "w") do file
        write(file, "t", t)
        write(file, "nx", result.nx)
        write(file, "ny", result.ny)
        write(file, "truncation_order", result.local_truncation_order)
        write(file, "rho_value", result.rho_value)
        write(file, "amplitude", result.amplitude)
        write(file, "dt", dt)
        write(file, "boundary_setup", Int(result.boundary_setup))

        integrator_names = [integrator_key_string(integ) for integ in result.integrators]
        write(file, "integrators", integrator_names)

        for integ in result.integrators
            key = integrator_key_string(integ)
            diag = result.diagnostics[integ]
            write(file, "density_mean_" * key, diag.rho_mean)
            write(file, "kinetic_energy_mean_" * key, diag.kinetic_energy_mean)
        end
    end

    return ts_path
end

function plot_lbe_integrator_comparison(result, output_path)
    integrators = result.integrators
    diagnostics = result.diagnostics
    t = dt .* collect(0:result.local_n_time - 1)
    baseline = integrators[1]
    baseline_diag = diagnostics[baseline]
    safelog(x) = max.(x, 1e-16)
    colors = ["k", "r", "b", "g", "m", "c"]

    close("all")
    fig = figure(figsize=(13, 10))

    subplot(2, 2, 1)
    for (idx, integ) in enumerate(integrators)
        diag = diagnostics[integ]
        plot(t, diag.rho_mean, color=colors[mod1(idx, length(colors))], linewidth=2.0, label=integrator_label(integ))
    end
    xlabel(L"t")
    ylabel(L"\langle \rho \rangle")
    title("Domain-averaged density")
    legend(loc="best")
    grid(true)

    subplot(2, 2, 2)
    for (idx, integ) in enumerate(integrators)
        diag = diagnostics[integ]
        semilogy(t, safelog(diag.kinetic_energy_mean), color=colors[mod1(idx, length(colors))], linewidth=2.0, label=integrator_label(integ))
    end
    xlabel(L"t")
    ylabel(L"\langle E_k \rangle")
    title("Domain-averaged kinetic energy")
    legend(loc="best")
    grid(true)

    subplot(2, 2, 3)
    for (idx, integ) in enumerate(integrators)
        diag = diagnostics[integ]
        semilogy(t, safelog(abs.(diag.rho_mean .- diag.rho_mean[1])), color=colors[mod1(idx, length(colors))], linewidth=2.0, label=integrator_label(integ))
    end
    xlabel(L"t")
    ylabel(L"|\langle \rho \rangle(t) - \langle \rho \rangle(0)|")
    title("Mean-density drift")
    legend(loc="best")
    grid(true)

    subplot(2, 2, 4)
    for (idx, integ) in enumerate(integrators)
        diag = diagnostics[integ]
        rel_to_baseline = abs.(diag.kinetic_energy_mean .- baseline_diag.kinetic_energy_mean)
        semilogy(t, safelog(rel_to_baseline), color=colors[mod1(idx, length(colors))], linewidth=2.0, label="$(integrator_label(integ)) vs $(integrator_label(baseline))")
    end
    xlabel(L"t")
    ylabel(L"|\langle E_k \rangle - \langle E_k \rangle_{\mathrm{baseline}}|")
    title("Kinetic-energy difference vs baseline")
    legend(loc="best")
    grid(true)

    boundary_label = result.boundary_setup ? "boundary-aware" : "periodic"
    suptitle("D2Q9 direct-LBE integrator comparison\n$(boundary_label), nx=$(result.nx), ny=$(result.ny), k=$(result.local_truncation_order), A=$(result.amplitude), ρ₀=$(result.rho_value), nt=$(result.local_n_time), dt=$(dt)")
    tight_layout(rect=(0, 0, 1, 0.95))
    savefig(output_path, bbox_inches="tight")
    println("Saved plot: $output_path")
end

function print_summary(result)
    println()
    println("Direct-LBE integrator comparison summary:")
    for integ in result.integrators
        diag = result.diagnostics[integ]
        @printf("  %-20s  mean(rho)[1]=%.8e  mean(rho)[end]=%.8e  mean(KE)[1]=%.8e  mean(KE)[end]=%.8e\n",
            integrator_label(integ),
            diag.rho_mean[1],
            diag.rho_mean[end],
            diag.kinetic_energy_mean[1],
            diag.kinetic_energy_mean[end],
        )
        @printf("  %-20s  max density drift = %.3e  max kinetic energy = %.3e\n",
            "",
            maximum(abs.(diag.rho_mean .- diag.rho_mean[1])),
            maximum(diag.kinetic_energy_mean),
        )
    end
end

function main(; nx=3,
    ny=3,
    amplitude=0.02,
    rho_value=1.0001,
    local_n_time=n_time,
    boundary_setup=false,
    coeff_method=:numerical,
    local_truncation_order=truncation_order,
    integrators=[:euler, :exponential_euler],
    save_output=true,
    output_dir=nothing,
    output_pdf=nothing)

    result = run_lbe_integrator_histories(
        nx=nx,
        ny=ny,
        amplitude=amplitude,
        rho_value=rho_value,
        local_n_time=local_n_time,
        boundary_setup=boundary_setup,
        coeff_method=coeff_method,
        local_truncation_order=local_truncation_order,
        integrators=integrators,
    )

    resolved_output_dir = isnothing(output_dir) ? default_output_dir(nx=nx, ny=ny, local_n_time=local_n_time, local_truncation_order=local_truncation_order) : output_dir
    mkpath(resolved_output_dir)
    resolved_output_pdf = isnothing(output_pdf) ? joinpath(resolved_output_dir, "d2q9_lbe_integrator_compare.pdf") : output_pdf

    plot_lbe_integrator_comparison(result, resolved_output_pdf)
    print_summary(result)

    saved_ts_path = nothing
    if save_output
        saved_ts_path = save_lbe_integrator_time_series_hdf5(result; output_dir=resolved_output_dir)
        println("Saved time series: $saved_ts_path")
    end

    return merge(result, (
        output_dir=resolved_output_dir,
        output_pdf=resolved_output_pdf,
        time_series_path=saved_ts_path,
    ))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end