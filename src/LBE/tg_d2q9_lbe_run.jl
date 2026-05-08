if !haskey(ENV, "QCFD_HOME")
    ENV["QCFD_HOME"] = normpath(joinpath(@__DIR__, "..", ".."))
end
if !haskey(ENV, "QCFD_SRC")
    ENV["QCFD_SRC"] = joinpath(ENV["QCFD_HOME"], "src") * "/"
end

QCFD_SRC = ENV["QCFD_SRC"]

using PyPlot
using LaTeXStrings
using LinearAlgebra

include(QCFD_SRC * "LBE/direct_LBE.jl")
include(QCFD_SRC * "LBE/d2q9_direct_setup.jl")

const DEFAULT_TEX_FIG_DIR = get(ENV, "QCFD_TEX_FIG_DIR", joinpath(homedir(), "Documents", "git-tex", "QC", "QCFD-QCLBM", "figs"))
const DEFAULT_TAU_VALUE = 1.0
const DEFAULT_DT = 0.1
const DEFAULT_N_TIME = 600
const DIRECT_LBE_POLY_ORDER = 3

function run_main_if_script(caller_file::AbstractString, main_fn::Function)
    if abspath(PROGRAM_FILE) == abspath(caller_file)
        main_fn()
    end
    return nothing
end

function field_l2_norm(ux, uy)
    return sqrt(sum(abs2, ux) + sum(abs2, uy))
end

function wall_speed_max(ux, uy)
    return maximum(abs.([ux[:, 1]; ux[:, end]; uy[:, 1]; uy[:, end]]))
end

function interior_speed_l2(ux, uy)
    if size(ux, 2) <= 2
        return 0.0
    end
    return sqrt(sum(abs2, ux[:, 2:end-1]) + sum(abs2, uy[:, 2:end-1]))
end

function speed_magnitude(ux, uy)
    return sqrt.(ux .^ 2 .+ uy .^ 2)
end

function default_plot_output_path(output_file, basename)
    if output_file !== nothing
        mkpath(dirname(output_file))
        return output_file
    end

    mkpath(DEFAULT_TEX_FIG_DIR)
    return joinpath(DEFAULT_TEX_FIG_DIR, basename * ".pdf")
end

function tg_decay_factor(nx, ny, tau_value, time_value)
    ν = (tau_value - 0.5) / 3
    kx = 2π / nx
    ky = 2π / ny
    return exp(-ν * (kx^2 + ky^2) * time_value)
end

function analytical_tg_velocity_field(i, j, nx, ny, amplitude, tau_value, time_value)
    ux0, uy0 = tg_velocity_field(i, j, nx, ny, amplitude)
    decay = tg_decay_factor(nx, ny, tau_value, time_value)
    return decay * ux0, decay * uy0
end

function analytical_tg_velocity_fields(nx, ny, amplitude, tau_value, time_value)
    ux = zeros(nx, ny)
    uy = zeros(nx, ny)
    for j in 1:ny
        for i in 1:nx
            ux[i, j], uy[i, j] = analytical_tg_velocity_field(i, j, nx, ny, amplitude, tau_value, time_value)
        end
    end
    return ux, uy
end

function velocity_error_metrics(ux_num, uy_num, ux_ref, uy_ref)
    err_vec = vcat(vec(ux_num - ux_ref), vec(uy_num - uy_ref))
    ref_vec = vcat(vec(ux_ref), vec(uy_ref))
    abs_l2 = norm(err_vec)
    rel_l2 = abs_l2 / max(norm(ref_vec), eps(Float64))
    max_abs = maximum(abs.([vec(ux_num - ux_ref); vec(uy_num - uy_ref)]))
    return abs_l2, rel_l2, max_abs
end

function plot_periodic_tg_analytical_comparison(ux_num, uy_num, ux_ref, uy_ref; output_file=nothing, tau_value, amplitude, final_time, abs_l2_hist=nothing, rel_l2_hist=nothing, max_abs_hist=nothing)
    speed_num = speed_magnitude(ux_num, uy_num)
    speed_ref = speed_magnitude(ux_ref, uy_ref)

    nx, ny = size(ux_num)
    mid_j = cld(ny, 2)
    x = collect(1:nx)
    have_error_histories = abs_l2_hist !== nothing && rel_l2_hist !== nothing && max_abs_hist !== nothing
    time_axis = have_error_histories ? collect(range(0.0, stop=final_time, length=length(abs_l2_hist))) : Float64[]
    X = repeat(reshape(collect(1:nx), :, 1), 1, ny)
    Y = repeat(reshape(collect(1:ny), 1, :), nx, 1)

    close("all")
    figure(figsize=have_error_histories ? (15, 8) : (15, 4.8))
    layout = have_error_histories ? (2, 3) : (1, 3)

    subplot(layout[1], layout[2], 1)
    quiver(X', Y', ux_num', uy_num')
    title("Numerical vectors")
    xlabel("x")
    ylabel("y")

    subplot(layout[1], layout[2], 2)
    imshow((speed_num - speed_ref)', origin="lower", cmap="RdBu")
    colorbar()
    title("|u| error")

    subplot(layout[1], layout[2], 3)
    plot(x, ux_num[:, mid_j], "or-", label="Direct LBE")
    plot(x, ux_ref[:, mid_j], "-k", linewidth=1.8, label="Analytical")
    xlabel("x")
    ylabel(L"u_x(x, y_{mid})")
    title("Centerline ux profile")
    legend(loc="best")

    if have_error_histories
        subplot(layout[1], layout[2], 4)
        plot(time_axis, abs_l2_hist, "-o", markersize=3, label="L2 abs")
        plot(time_axis, max_abs_hist, "-s", markersize=3, label="Max abs")
        xlabel("t")
        ylabel("Error")
        title("Absolute error vs time")
        legend(loc="best")

        subplot(layout[1], layout[2], 5)
        plot(time_axis, rel_l2_hist, "-o", markersize=3, color="tab:red")
        xlabel("t")
        ylabel("Relative L2 error")
        title("Relative error vs time")
    end

    suptitle("Periodic D2Q9 TG direct LBE vs analytical, τ=$(tau_value), amplitude=$(amplitude), t=$(final_time)")
    tight_layout(rect=(0, 0, 1, 0.96))

    if output_file !== nothing
        savefig(output_file, dpi=200, bbox_inches="tight")
    end

    display(gcf())
    show()
end

function plot_tg_boundary_velocity_fields(ux_initial, uy_initial, ux_final, uy_final; output_file=nothing, tau_value, amplitude, n_time)
    speed_initial = speed_magnitude(ux_initial, uy_initial)
    speed_final = speed_magnitude(ux_final, uy_final)

    nx, ny = size(ux_initial)
    x = collect(1:nx)
    y = collect(1:ny)
    X = repeat(reshape(x, :, 1), 1, ny)
    Y = repeat(reshape(y, 1, :), nx, 1)

    close("all")
    figure(figsize=(14, 8))

    subplot(2, 3, 1)
    imshow(speed_initial', origin="lower", cmap="viridis")
    colorbar()
    title("Initial |u|")

    subplot(2, 3, 2)
    quiver(X', Y', ux_initial', uy_initial')
    title("Initial velocity vectors")
    xlabel("x")
    ylabel("y")

    subplot(2, 3, 3)
    imshow(ux_initial', origin="lower", cmap="RdBu")
    colorbar()
    title(L"Initial\ u_x")

    subplot(2, 3, 4)
    imshow(speed_final', origin="lower", cmap="viridis")
    colorbar()
    title("Final |u|")

    subplot(2, 3, 5)
    quiver(X', Y', ux_final', uy_final')
    title("Final velocity vectors")
    xlabel("x")
    ylabel("y")

    subplot(2, 3, 6)
    imshow(ux_final', origin="lower", cmap="RdBu")
    colorbar()
    title(L"Final\ u_x")

    suptitle("D2Q9 TG boundary-value direct LBE, τ=$(tau_value), amplitude=$(amplitude), n_t=$(n_time)")
    tight_layout(rect=(0, 0, 1, 0.96))

    if output_file !== nothing
        savefig(output_file, dpi=200, bbox_inches="tight")
    end

    display(gcf())
    show()
end

function macroscopic_histories_from_phi_history(phi_hist, nx, ny, e_value)
    n_time_local = size(phi_hist, 2)
    rho_hist = zeros(nx, ny, n_time_local)
    ux_hist = zeros(nx, ny, n_time_local)
    uy_hist = zeros(nx, ny, n_time_local)

    for nt in 1:n_time_local
        rho, ux, uy = macroscopic_fields_from_state(phi_hist[:, nt], nx, ny, e_value)
        rho_hist[:, :, nt] .= rho
        ux_hist[:, :, nt] .= ux
        uy_hist[:, :, nt] .= uy
    end

    return rho_hist, ux_hist, uy_hist
end

function analytical_velocity_error_histories_dt(ux_hist, uy_hist, nx, ny, amplitude, tau_value_local, dt_local)
    n_time_local = size(ux_hist, 3)
    abs_l2_hist = zeros(n_time_local)
    rel_l2_hist = zeros(n_time_local)
    max_abs_hist = zeros(n_time_local)

    for nt in 1:n_time_local
        time_value = dt_local * (nt - 1)
        ux_ref, uy_ref = analytical_tg_velocity_fields(nx, ny, amplitude, tau_value_local, time_value)
        abs_l2_hist[nt], rel_l2_hist[nt], max_abs_hist[nt] = velocity_error_metrics(
            ux_hist[:, :, nt],
            uy_hist[:, :, nt],
            ux_ref,
            uy_ref,
        )
    end

    return abs_l2_hist, rel_l2_hist, max_abs_hist
end

function direct_lbe_tg_output_file(output_file, basename)
    return default_plot_output_path(output_file, basename)
end

function run_tg_d2q9_lbe_core(; nx=3,
    ny=3,
    amplitude=0.02,
    rho_value=1.0001,
    tau_value_input=DEFAULT_TAU_VALUE,
    dt_value=DEFAULT_DT,
    local_n_time=DEFAULT_N_TIME,
    boundary_setup=false,
    direct_lbe_integrator=:euler,
    return_phi_history=false,
    l_plot=false,
    output_file=nothing,
    compare_analytical=false)

    runtime = build_direct_lbe_d2q9_runtime(
        nx=nx,
        ny=ny,
        rho_value=rho_value,
        tau_value=tau_value_input,
        boundary_setup=boundary_setup,
        poly_order_input=DIRECT_LBE_POLY_ORDER,
        lTaylor_input=true,
    )

    reference = build_direct_lbe_tg_reference(
        nx=nx,
        ny=ny,
        amplitude=amplitude,
        rho_value=rho_value,
        local_n_time=local_n_time,
        boundary_setup=boundary_setup,
        runtime=runtime,
        dt_value=dt_value,
        direct_lbe_integrator=direct_lbe_integrator,
    )

    phi_hist = reference.reference_phi_history
    rho_hist, ux_hist, uy_hist = macroscopic_histories_from_phi_history(phi_hist, nx, ny, runtime.e_value)
    mass_hist = vec([sum(rho_hist[:, :, nt]) for nt in 1:local_n_time])
    speed_norm_hist = vec([field_l2_norm(ux_hist[:, :, nt], uy_hist[:, :, nt]) for nt in 1:local_n_time])
    wall_speed_hist = vec([wall_speed_max(ux_hist[:, :, nt], uy_hist[:, :, nt]) for nt in 1:local_n_time])
    interior_speed_hist = vec([interior_speed_l2(ux_hist[:, :, nt], uy_hist[:, :, nt]) for nt in 1:local_n_time])

    println(boundary_setup ? "Running boundary-aware D2Q9 TG direct LBE" : "Running periodic D2Q9 TG direct LBE")
    println("  grid = $(nx)×$(ny)")
    println("  tau = $tau_value_input")
    println("  dt = $dt_value")
    println("  amplitude = $amplitude")
    println("  rho_value = $rho_value")
    println("  direct_lbe_integrator = $(normalize_direct_lbe_integrator(direct_lbe_integrator))")
    println("  initial total mass = ", mass_hist[1])
    println("  final total mass = ", mass_hist[end])
    println("  max |mass(t) - mass(0)| = ", maximum(abs.(mass_hist .- mass_hist[1])))
    println("  initial velocity L2 norm = ", speed_norm_hist[1])
    println("  final velocity L2 norm = ", speed_norm_hist[end])
    if boundary_setup
        println("  initial wall max speed = ", wall_speed_hist[1])
        println("  final wall max speed = ", wall_speed_hist[end])
        println("  initial interior velocity L2 norm = ", interior_speed_hist[1])
        println("  final interior velocity L2 norm = ", interior_speed_hist[end])
    end
    println("  initial ux field = ", ux_hist[:, :, 1])
    println("  initial uy field = ", uy_hist[:, :, 1])
    println("  final ux field = ", ux_hist[:, :, end])
    println("  final uy field = ", uy_hist[:, :, end])

    resolved_output_file = nothing
    if l_plot
        if boundary_setup
            resolved_output_file = direct_lbe_tg_output_file(output_file, "tg_d2q9_lbe_boundary_nx$(nx)_ny$(ny)_nt$(local_n_time)")
        elseif compare_analytical
            resolved_output_file = direct_lbe_tg_output_file(output_file, "tg_d2q9_lbe_periodic_vs_analytical_nx$(nx)_ny$(ny)_nt$(local_n_time)")
        else
            resolved_output_file = direct_lbe_tg_output_file(output_file, "tg_d2q9_lbe_periodic_nx$(nx)_ny$(ny)_nt$(local_n_time)")
        end
    end

    abs_l2_hist = nothing
    rel_l2_hist = nothing
    max_abs_hist = nothing
    if compare_analytical && !boundary_setup
        abs_l2_hist, rel_l2_hist, max_abs_hist = analytical_velocity_error_histories_dt(
            ux_hist,
            uy_hist,
            nx,
            ny,
            amplitude,
            tau_value_input,
            dt_value,
        )
        final_time = dt_value * (local_n_time - 1)
        ux_ref, uy_ref = analytical_tg_velocity_fields(nx, ny, amplitude, tau_value_input, final_time)
        abs_l2, rel_l2, max_abs = velocity_error_metrics(ux_hist[:, :, end], uy_hist[:, :, end], ux_ref, uy_ref)
        println("  analytical comparison time = ", final_time)
        println("  analytical velocity L2 error = ", abs_l2)
        println("  analytical velocity relative L2 error = ", rel_l2)
        println("  analytical velocity max abs error = ", max_abs)

        if l_plot
            plot_periodic_tg_analytical_comparison(
                ux_hist[:, :, end],
                uy_hist[:, :, end],
                ux_ref,
                uy_ref;
                output_file=resolved_output_file,
                tau_value=tau_value_input,
                amplitude=amplitude,
                final_time=final_time,
                abs_l2_hist=abs_l2_hist,
                rel_l2_hist=rel_l2_hist,
                max_abs_hist=max_abs_hist,
            )
        end
    elseif l_plot && boundary_setup
        plot_tg_boundary_velocity_fields(
            ux_hist[:, :, 1],
            uy_hist[:, :, 1],
            ux_hist[:, :, end],
            uy_hist[:, :, end];
            output_file=resolved_output_file,
            tau_value=tau_value_input,
            amplitude=amplitude,
            n_time=local_n_time,
        )
    end

    if resolved_output_file !== nothing
        println("  saved figure = ", resolved_output_file)
    end

    result = (
        rho_hist=rho_hist,
        ux_hist=ux_hist,
        uy_hist=uy_hist,
        mass_hist=mass_hist,
        speed_norm_hist=speed_norm_hist,
        wall_speed_hist=wall_speed_hist,
        interior_speed_hist=interior_speed_hist,
        phi_hist=phi_hist,
        e_value=runtime.e_value,
        output_file=resolved_output_file,
        analytical_abs_l2_hist=abs_l2_hist,
        analytical_rel_l2_hist=rel_l2_hist,
        analytical_max_abs_hist=max_abs_hist,
    )

    if return_phi_history
        return result
    end

    return merge(result, (phi_hist=nothing,))
end

function run_tg_d2q9_lbe(; nx=3,
    ny=3,
    amplitude=0.02,
    rho_value=1.0001,
    tau_value=DEFAULT_TAU_VALUE,
    dt=DEFAULT_DT,
    n_time=DEFAULT_N_TIME,
    return_phi_history=false,
    l_plot=false,
    output_file=nothing,
    compare_analytical=false,
    direct_lbe_integrator=:euler)

    result = run_tg_d2q9_lbe_core(
        nx=nx,
        ny=ny,
        amplitude=amplitude,
        rho_value=rho_value,
        tau_value_input=tau_value,
        dt_value=dt,
        local_n_time=n_time,
        boundary_setup=false,
        direct_lbe_integrator=direct_lbe_integrator,
        return_phi_history=return_phi_history,
        l_plot=l_plot,
        output_file=output_file,
        compare_analytical=compare_analytical,
    )

    if return_phi_history
        return result.rho_hist, result.ux_hist, result.uy_hist, result.mass_hist, result.speed_norm_hist, result.phi_hist
    end
    return result.rho_hist, result.ux_hist, result.uy_hist, result.mass_hist, result.speed_norm_hist
end

function run_tg_d2q9_boundary_lbe(; nx=3,
    ny=3,
    amplitude=0.02,
    rho_value=1.0001,
    tau_value=DEFAULT_TAU_VALUE,
    dt=DEFAULT_DT,
    n_time=DEFAULT_N_TIME,
    return_phi_history=false,
    l_plot=false,
    output_file=nothing,
    direct_lbe_integrator=:euler)

    result = run_tg_d2q9_lbe_core(
        nx=nx,
        ny=ny,
        amplitude=amplitude,
        rho_value=rho_value,
        tau_value_input=tau_value,
        dt_value=dt,
        local_n_time=n_time,
        boundary_setup=true,
        direct_lbe_integrator=direct_lbe_integrator,
        return_phi_history=return_phi_history,
        l_plot=l_plot,
        output_file=output_file,
        compare_analytical=false,
    )

    if return_phi_history
        return result.rho_hist, result.ux_hist, result.uy_hist, result.mass_hist, result.speed_norm_hist, result.wall_speed_hist, result.interior_speed_hist, result.phi_hist
    end
    return result.rho_hist, result.ux_hist, result.uy_hist, result.mass_hist, result.speed_norm_hist, result.wall_speed_hist, result.interior_speed_hist
end

function main(; nx=8,
    ny=8,
    amplitude=0.02,
    rho_value=1.0001,
    tau_value=DEFAULT_TAU_VALUE,
    dt=DEFAULT_DT,
    n_time=50,
    boundary_setup=false,
    direct_lbe_integrator=:euler,
    compare_analytical=!boundary_setup,
    l_plot=true,
    output_file=nothing)

    if boundary_setup
        return run_tg_d2q9_boundary_lbe(
            nx=nx,
            ny=ny,
            amplitude=amplitude,
            rho_value=rho_value,
            tau_value=tau_value,
            dt=dt,
            n_time=n_time,
            l_plot=l_plot,
            output_file=output_file,
            direct_lbe_integrator=direct_lbe_integrator,
        )
    end

    return run_tg_d2q9_lbe(
        nx=nx,
        ny=ny,
        amplitude=amplitude,
        rho_value=rho_value,
        tau_value=tau_value,
        dt=dt,
        n_time=n_time,
        l_plot=l_plot,
        output_file=output_file,
        compare_analytical=compare_analytical,
        direct_lbe_integrator=direct_lbe_integrator,
    )
end

run_main_if_script(@__FILE__, main)