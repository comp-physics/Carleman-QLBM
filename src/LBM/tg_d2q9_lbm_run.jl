QCFD_SRC = ENV["QCFD_SRC"]

using PyPlot
using LaTeXStrings
using LinearAlgebra

include(QCFD_SRC * "LBM/lbm_cons.jl")
include(QCFD_SRC * "LBM/cal_feq.jl")
include(QCFD_SRC * "LBM/streaming.jl")

const DEFAULT_TEX_FIG_DIR = get(ENV, "QCFD_TEX_FIG_DIR", joinpath(homedir(), "Documents", "git-tex", "QC", "QCFD-QCLBM", "figs"))

function tg_velocity_field_theory(i, j, nx, ny, amplitude)
    rx = -π + 2π * (i - 1) / nx
    ry = -π + 2π * (j - 1) / ny
    ux = amplitude * sin(rx) * cos(ry)
    uy = -amplitude * cos(rx) * sin(ry)
    return ux, uy
end

function equilibrium_distribution(rho, ux, uy, weights, velocities, a, b, c, d)
    vcx = velocities[:, 1]
    vcy = velocities[:, 2]
    return cal_feq(rho, ux, uy, weights, vcx, vcy, a, b, c, d)
end

function flatten_lbd_state(LBD_state)
    nx, ny, q = size(LBD_state)
    phi = zeros(nx * ny * q)

    for j in 1:ny
        for i in 1:nx
            start_idx = ((j - 1) * nx + (i - 1)) * q + 1
            phi[start_idx:start_idx + q - 1] .= vec(LBD_state[i, j, :])
        end
    end

    return phi
end

function shifted_incompressible_equilibrium(delta_rho, ux, uy, weights, velocities)
    g_eq = zeros(length(weights))
    u2 = ux^2 + uy^2

    for m in eachindex(weights)
        eu = velocities[m, 1] * ux + velocities[m, 2] * uy
        g_eq[m] = weights[m] * (delta_rho + 3 * eu + 4.5 * eu^2 - 1.5 * u2)
    end

    return g_eq
end

function macroscopic_fields_from_lbd(LBD_state, velocities)
    nx, ny, q = size(LBD_state)
    rho = zeros(nx, ny)
    ux = zeros(nx, ny)
    uy = zeros(nx, ny)

    for j in 1:ny
        for i in 1:nx
            fi = vec(LBD_state[i, j, :])
            rho_local = sum(fi)
            rho[i, j] = rho_local
            ux[i, j] = sum(fi[k] * velocities[k, 1] for k = 1:q) / rho_local
            uy[i, j] = sum(fi[k] * velocities[k, 2] for k = 1:q) / rho_local
        end
    end

    return rho, ux, uy
end

function initialize_tg_lbd(nx, ny, amplitude, rho_value, weights, velocities, a, b, c, d)
    q = length(weights)
    LBD = zeros(nx, ny, q, 2)

    for j in 1:ny
        for i in 1:nx
            ux, uy = tg_velocity_field_theory(i, j, nx, ny, amplitude)
            feq = equilibrium_distribution(rho_value, ux, uy, weights, velocities, a, b, c, d)
            LBD[i, j, :, 1] .= feq
            LBD[i, j, :, 2] .= feq
        end
    end

    return LBD
end

function initialize_tg_boundary_lbd(nx, ny, amplitude, rho_value, weights, velocities, a, b, c, d)
    q = length(weights)
    LBD = zeros(nx, ny, q, 2)

    for j in 1:ny
        for i in 1:nx
            if j == 1 || j == ny
                ux = 0.0
                uy = 0.0
            else
                ux, uy = tg_velocity_field_theory(i, j, nx, ny, amplitude)
            end

            feq = equilibrium_distribution(rho_value, ux, uy, weights, velocities, a, b, c, d)
            LBD[i, j, :, 1] .= feq
            LBD[i, j, :, 2] .= feq
        end
    end

    return LBD
end

function collide!(LBD, tau_value, weights, velocities, a, b, c, d)
    nx, ny, q = size(LBD[:, :, :, 1])
    rho, ux, uy = macroscopic_fields_from_lbd(LBD[:, :, :, 1], velocities)

    for j in 1:ny
        for i in 1:nx
            feq = equilibrium_distribution(rho[i, j], ux[i, j], uy[i, j], weights, velocities, a, b, c, d)
            for k in 1:q
                LBD[i, j, k, 2] = LBD[i, j, k, 1] - (LBD[i, j, k, 1] - feq[k]) / tau_value
            end
        end
    end

    return LBD
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
    ν = lbm_param() * (tau_value - 0.5)
    kx = 2π / nx
    ky = 2π / ny
    return exp(-ν * (kx^2 + ky^2) * time_value)
end

function analytical_tg_velocity_field(i, j, nx, ny, amplitude, tau_value, time_value)
    ux0, uy0 = tg_velocity_field_theory(i, j, nx, ny, amplitude)
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

function analytical_velocity_error_histories(ux_hist, uy_hist, nx, ny, amplitude, tau_value)
    n_time = size(ux_hist, 3)
    abs_l2_hist = zeros(n_time)
    rel_l2_hist = zeros(n_time)
    max_abs_hist = zeros(n_time)

    for nt in 1:n_time
        ux_ref, uy_ref = analytical_tg_velocity_fields(nx, ny, amplitude, tau_value, nt - 1)
        abs_l2_hist[nt], rel_l2_hist[nt], max_abs_hist[nt] = velocity_error_metrics(
            ux_hist[:, :, nt],
            uy_hist[:, :, nt],
            ux_ref,
            uy_ref,
        )
    end

    return abs_l2_hist, rel_l2_hist, max_abs_hist
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
    plot(x, ux_num[:, mid_j], "or-", label="LBM")
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

    suptitle("Periodic D2Q9 TG vs analytical, τ=$(tau_value), amplitude=$(amplitude), t=$(final_time)")
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

    suptitle("D2Q9 TG boundary-value LBM, τ=$(tau_value), amplitude=$(amplitude), n_t=$(n_time)")
    tight_layout(rect=(0, 0, 1, 0.96))

    if output_file !== nothing
        savefig(output_file, dpi=200, bbox_inches="tight")
    end

    display(gcf())
    show()
end

function macroscopic_fields_from_shifted_lbd(LBD_state, velocities)
    nx, ny, q = size(LBD_state)
    delta_rho = zeros(nx, ny)
    ux = zeros(nx, ny)
    uy = zeros(nx, ny)

    for j in 1:ny
        for i in 1:nx
            gi = vec(LBD_state[i, j, :])
            delta_rho_local = sum(gi)
            delta_rho[i, j] = delta_rho_local
            ux[i, j] = sum(gi[k] * velocities[k, 1] for k = 1:q)
            uy[i, j] = sum(gi[k] * velocities[k, 2] for k = 1:q)
        end
    end

    return delta_rho, ux, uy
end

function initialize_tg_shifted_lbd(nx, ny, amplitude, weights, velocities)
    q = length(weights)
    LBD = zeros(nx, ny, q, 2)

    for j in 1:ny
        for i in 1:nx
            ux, uy = tg_velocity_field_theory(i, j, nx, ny, amplitude)
            g_eq = shifted_incompressible_equilibrium(0.0, ux, uy, weights, velocities)
            LBD[i, j, :, 1] .= g_eq
            LBD[i, j, :, 2] .= g_eq
        end
    end

    return LBD
end

function collide_shifted!(LBD, tau_value, weights, velocities)
    nx, ny, q = size(LBD[:, :, :, 1])
    delta_rho, ux, uy = macroscopic_fields_from_shifted_lbd(LBD[:, :, :, 1], velocities)

    for j in 1:ny
        for i in 1:nx
            g_eq = shifted_incompressible_equilibrium(delta_rho[i, j], ux[i, j], uy[i, j], weights, velocities)
            for k in 1:q
                LBD[i, j, k, 2] = LBD[i, j, k, 1] - (LBD[i, j, k, 1] - g_eq[k]) / tau_value
            end
        end
    end

    return LBD
end

function run_tg_d2q9_lbm(; nx=3, ny=3, amplitude=0.02, rho_value=1.0, tau_value=0.8, n_time=5, l_noslipBC=false, return_phi_history=false, l_plot=false, output_file=nothing, compare_analytical=false)
    global LX = nx
    global LY = ny
    a, b, c, d, D_local, Q_local, weights, velocities = lbm_cons()
    global D = D_local
    global Q = Q_local
    global e = velocities

    LBD = initialize_tg_lbd(nx, ny, amplitude, rho_value, weights, velocities, a, b, c, d)

    rho_hist = zeros(nx, ny, n_time)
    ux_hist = zeros(nx, ny, n_time)
    uy_hist = zeros(nx, ny, n_time)
    mass_hist = zeros(n_time)
    speed_norm_hist = zeros(n_time)
    phi_hist = return_phi_history ? zeros(nx * ny * Q_local, n_time) : nothing

    rho0, ux0, uy0 = macroscopic_fields_from_lbd(LBD[:, :, :, 1], velocities)
    rho_hist[:, :, 1] .= rho0
    ux_hist[:, :, 1] .= ux0
    uy_hist[:, :, 1] .= uy0
    mass_hist[1] = sum(rho0)
    speed_norm_hist[1] = field_l2_norm(ux0, uy0)
    if return_phi_history
        phi_hist[:, 1] .= flatten_lbd_state(LBD[:, :, :, 1])
    end

    for nt in 2:n_time
        collide!(LBD, tau_value, weights, velocities, a, b, c, d)
        LBD, _ = streaming(LBD, nx, ny, l_noslipBC)

        rho, ux, uy = macroscopic_fields_from_lbd(LBD[:, :, :, 1], velocities)
        rho_hist[:, :, nt] .= rho
        ux_hist[:, :, nt] .= ux
        uy_hist[:, :, nt] .= uy
        mass_hist[nt] = sum(rho)
        speed_norm_hist[nt] = field_l2_norm(ux, uy)
        if return_phi_history
            phi_hist[:, nt] .= flatten_lbd_state(LBD[:, :, :, 1])
        end
    end

    println("Running pure D2Q9 TG LBM")
    println("  grid = $(nx)×$(ny)")
    println("  tau = $tau_value")
    println("  amplitude = $amplitude")
    println("  no-slip BC = $l_noslipBC")
    println("  initial total mass = ", mass_hist[1])
    println("  final total mass = ", mass_hist[end])
    println("  max |mass(t) - mass(0)| = ", maximum(abs.(mass_hist .- mass_hist[1])))
    println("  initial velocity L2 norm = ", speed_norm_hist[1])
    println("  final velocity L2 norm = ", speed_norm_hist[end])
    println("  initial ux field = ", ux_hist[:, :, 1])
    println("  initial uy field = ", uy_hist[:, :, 1])
    println("  final ux field = ", ux_hist[:, :, end])
    println("  final uy field = ", uy_hist[:, :, end])

    resolved_output_file = nothing
    if l_plot
        if compare_analytical
            resolved_output_file = default_plot_output_path(
                output_file,
                "tg_d2q9_periodic_vs_analytical_nx$(nx)_ny$(ny)_nt$(n_time)",
            )
        else
            resolved_output_file = default_plot_output_path(
                output_file,
                "tg_d2q9_periodic_nx$(nx)_ny$(ny)_nt$(n_time)",
            )
        end
    end

    if compare_analytical
        abs_l2_hist, rel_l2_hist, max_abs_hist = analytical_velocity_error_histories(
            ux_hist,
            uy_hist,
            nx,
            ny,
            amplitude,
            tau_value,
        )
        final_time = n_time - 1
        ux_ref, uy_ref = analytical_tg_velocity_fields(nx, ny, amplitude, tau_value, final_time)
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
                tau_value=tau_value,
                amplitude=amplitude,
                final_time=final_time,
                abs_l2_hist=abs_l2_hist,
                rel_l2_hist=rel_l2_hist,
                max_abs_hist=max_abs_hist,
            )
        end
    end

    if resolved_output_file !== nothing
        println("  saved figure = ", resolved_output_file)
    end

    if return_phi_history
        return rho_hist, ux_hist, uy_hist, mass_hist, speed_norm_hist, phi_hist
    end

    return rho_hist, ux_hist, uy_hist, mass_hist, speed_norm_hist
end

function run_tg_d2q9_boundary_lbm(; nx=3, ny=3, amplitude=0.02, rho_value=1.0, tau_value=0.8, n_time=5, return_phi_history=false, l_plot=false, output_file=nothing)
    global LX = nx
    global LY = ny
    a, b, c, d, D_local, Q_local, weights, velocities = lbm_cons()
    global D = D_local
    global Q = Q_local
    global e = velocities

    LBD = initialize_tg_boundary_lbd(nx, ny, amplitude, rho_value, weights, velocities, a, b, c, d)

    rho_hist = zeros(nx, ny, n_time)
    ux_hist = zeros(nx, ny, n_time)
    uy_hist = zeros(nx, ny, n_time)
    mass_hist = zeros(n_time)
    speed_norm_hist = zeros(n_time)
    wall_speed_hist = zeros(n_time)
    interior_speed_hist = zeros(n_time)
    phi_hist = return_phi_history ? zeros(nx * ny * Q_local, n_time) : nothing

    rho0, ux0, uy0 = macroscopic_fields_from_lbd(LBD[:, :, :, 1], velocities)
    rho_hist[:, :, 1] .= rho0
    ux_hist[:, :, 1] .= ux0
    uy_hist[:, :, 1] .= uy0
    mass_hist[1] = sum(rho0)
    speed_norm_hist[1] = field_l2_norm(ux0, uy0)
    wall_speed_hist[1] = wall_speed_max(ux0, uy0)
    interior_speed_hist[1] = interior_speed_l2(ux0, uy0)
    if return_phi_history
        phi_hist[:, 1] .= flatten_lbd_state(LBD[:, :, :, 1])
    end

    for nt in 2:n_time
        collide!(LBD, tau_value, weights, velocities, a, b, c, d)
        LBD, _ = streaming(LBD, nx, ny, true)

        rho, ux, uy = macroscopic_fields_from_lbd(LBD[:, :, :, 1], velocities)
        rho_hist[:, :, nt] .= rho
        ux_hist[:, :, nt] .= ux
        uy_hist[:, :, nt] .= uy
        mass_hist[nt] = sum(rho)
        speed_norm_hist[nt] = field_l2_norm(ux, uy)
        wall_speed_hist[nt] = wall_speed_max(ux, uy)
        interior_speed_hist[nt] = interior_speed_l2(ux, uy)
        if return_phi_history
            phi_hist[:, nt] .= flatten_lbd_state(LBD[:, :, :, 1])
        end
    end

    println("Running boundary-value D2Q9 TG LBM with current scheme")
    println("  grid = $(nx)×$(ny)")
    println("  tau = $tau_value")
    println("  amplitude = $amplitude")
    println("  top/bottom no-slip walls = true")
    println("  initial total mass = ", mass_hist[1])
    println("  final total mass = ", mass_hist[end])
    println("  max |mass(t) - mass(0)| = ", maximum(abs.(mass_hist .- mass_hist[1])))
    println("  initial full velocity L2 norm = ", speed_norm_hist[1])
    println("  final full velocity L2 norm = ", speed_norm_hist[end])
    println("  initial wall max speed = ", wall_speed_hist[1])
    println("  final wall max speed = ", wall_speed_hist[end])
    println("  initial interior velocity L2 norm = ", interior_speed_hist[1])
    println("  final interior velocity L2 norm = ", interior_speed_hist[end])
    println("  initial ux field = ", ux_hist[:, :, 1])
    println("  initial uy field = ", uy_hist[:, :, 1])
    println("  final ux field = ", ux_hist[:, :, end])
    println("  final uy field = ", uy_hist[:, :, end])

    resolved_output_file = l_plot ? default_plot_output_path(
        output_file,
        "tg_d2q9_boundary_nx$(nx)_ny$(ny)_nt$(n_time)",
    ) : nothing

    if l_plot
        plot_tg_boundary_velocity_fields(
            ux_hist[:, :, 1],
            uy_hist[:, :, 1],
            ux_hist[:, :, end],
            uy_hist[:, :, end];
            output_file=resolved_output_file,
            tau_value=tau_value,
            amplitude=amplitude,
            n_time=n_time,
        )
    end

    if resolved_output_file !== nothing
        println("  saved figure = ", resolved_output_file)
    end

    if return_phi_history
        return rho_hist, ux_hist, uy_hist, mass_hist, speed_norm_hist, wall_speed_hist, interior_speed_hist, phi_hist
    end

    return rho_hist, ux_hist, uy_hist, mass_hist, speed_norm_hist, wall_speed_hist, interior_speed_hist
end

function run_tg_d2q9_shifted_lbm(; nx=3, ny=3, amplitude=0.02, tau_value=0.8, n_time=5, l_noslipBC=false)
    global LX = nx
    global LY = ny
    a, b, c, d, D_local, Q_local, weights, velocities = lbm_cons()
    global D = D_local
    global Q = Q_local
    global e = velocities

    LBD = initialize_tg_shifted_lbd(nx, ny, amplitude, weights, velocities)

    delta_rho_hist = zeros(nx, ny, n_time)
    ux_hist = zeros(nx, ny, n_time)
    uy_hist = zeros(nx, ny, n_time)
    delta_mass_hist = zeros(n_time)
    speed_norm_hist = zeros(n_time)

    delta_rho0, ux0, uy0 = macroscopic_fields_from_shifted_lbd(LBD[:, :, :, 1], velocities)
    delta_rho_hist[:, :, 1] .= delta_rho0
    ux_hist[:, :, 1] .= ux0
    uy_hist[:, :, 1] .= uy0
    delta_mass_hist[1] = sum(delta_rho0)
    speed_norm_hist[1] = field_l2_norm(ux0, uy0)

    for nt in 2:n_time
        collide_shifted!(LBD, tau_value, weights, velocities)
        LBD, _ = streaming(LBD, nx, ny, l_noslipBC)

        delta_rho, ux, uy = macroscopic_fields_from_shifted_lbd(LBD[:, :, :, 1], velocities)
        delta_rho_hist[:, :, nt] .= delta_rho
        ux_hist[:, :, nt] .= ux
        uy_hist[:, :, nt] .= uy
        delta_mass_hist[nt] = sum(delta_rho)
        speed_norm_hist[nt] = field_l2_norm(ux, uy)
    end

    println("Running shifted-incompressible D2Q9 TG LBM")
    println("  grid = $(nx)×$(ny)")
    println("  tau = $tau_value")
    println("  amplitude = $amplitude")
    println("  no-slip BC = $l_noslipBC")
    println("  initial total δρ = ", delta_mass_hist[1])
    println("  final total δρ = ", delta_mass_hist[end])
    println("  max |Σδρ(t) - Σδρ(0)| = ", maximum(abs.(delta_mass_hist .- delta_mass_hist[1])))
    println("  max |δρ| = ", maximum(abs.(delta_rho_hist)))
    println("  initial velocity L2 norm = ", speed_norm_hist[1])
    println("  final velocity L2 norm = ", speed_norm_hist[end])
    println("  initial ux field = ", ux_hist[:, :, 1])
    println("  initial uy field = ", uy_hist[:, :, 1])
    println("  final ux field = ", ux_hist[:, :, end])
    println("  final uy field = ", uy_hist[:, :, end])

    return delta_rho_hist, ux_hist, uy_hist, delta_mass_hist, speed_norm_hist
end