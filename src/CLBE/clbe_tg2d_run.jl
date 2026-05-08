l_sympy = true
QCFD_SRC = ENV["QCFD_SRC"]
QCFD_HOME = ENV["QCFD_HOME"]

using PyPlot
using LaTeXStrings
using HDF5
using Statistics
using Printf

include(QCFD_HOME * "/visualization/plot_kit.jl")

if l_sympy
    using SymPy
    using LinearAlgebra
    using SparseArrays
    include(QCFD_SRC * "CLBE/coeffs_poly.jl")
else
    using Symbolics
end

include("clbe_config_2D.jl")

include(QCFD_SRC * "CLBE/collision_sym.jl")
include(QCFD_SRC * "CLBE/carleman_transferA.jl")
include(QCFD_SRC * "CLBE/carleman_transferA_ngrid.jl")
include(QCFD_SRC * "CLBE/LBM_const_subs.jl")
include(QCFD_SRC * "LBM/lbm_const_sym.jl")
include(QCFD_SRC * "LBM/cal_feq.jl")
include(QCFD_SRC * "LBM/tg_d2q9_lbm_run.jl")
include(QCFD_SRC * "LBE/direct_LBE.jl")
include(QCFD_SRC * "CLBE/timeMarching.jl")

function streaming_operator_D2Q9_interleaved_periodic(nx, ny, hx, hy)
    e = [
        [0, 0],
        [1, 0],
        [-1, 0],
        [0, 1],
        [0, -1],
        [1, 1],
        [1, -1],
        [-1, 1],
        [-1, -1],
    ]

    n_velocities = 9
    n_total = n_velocities * nx * ny

    I_idx = Int[]
    J_idx = Int[]
    V_vals = Float64[]

    global_index(vel_idx, i, j) = ((j - 1) * nx + (i - 1)) * n_velocities + vel_idx

    for j in 1:ny
        for i in 1:nx
            for vel in 1:n_velocities
                row_idx = global_index(vel, i, j)
                ex, ey = e[vel]

                if ex == 0 && ey == 0
                    push!(I_idx, row_idx)
                    push!(J_idx, row_idx)
                    push!(V_vals, 0.0)
                    continue
                end

                if ex != 0
                    left_i = i == 1 ? nx : i - 1
                    right_i = i == nx ? 1 : i + 1

                    push!(I_idx, row_idx)
                    push!(J_idx, global_index(vel, left_i, j))
                    push!(V_vals, -ex / (2 * hx))

                    push!(I_idx, row_idx)
                    push!(J_idx, global_index(vel, right_i, j))
                    push!(V_vals, ex / (2 * hx))
                end

                if ey != 0
                    bottom_j = j == 1 ? ny : j - 1
                    top_j = j == ny ? 1 : j + 1

                    push!(I_idx, row_idx)
                    push!(J_idx, global_index(vel, i, bottom_j))
                    push!(V_vals, -ey / (2 * hy))

                    push!(I_idx, row_idx)
                    push!(J_idx, global_index(vel, i, top_j))
                    push!(V_vals, ey / (2 * hy))
                end
            end
        end
    end

    return sparse(I_idx, J_idx, V_vals, n_total, n_total), e
end

function opposite_velocity_index_D2Q9(vel)
    opposite_map = Dict(
        1 => 1,
        2 => 3,
        3 => 2,
        4 => 5,
        5 => 4,
        6 => 9,
        7 => 8,
        8 => 7,
        9 => 6,
    )
    return opposite_map[vel]
end

function streaming_operator_D2Q9_interleaved_boundary_aware(nx, ny, hx, hy)
    e = [
        [0, 0],
        [1, 0],
        [-1, 0],
        [0, 1],
        [0, -1],
        [1, 1],
        [1, -1],
        [-1, 1],
        [-1, -1],
    ]

    n_velocities = 9
    n_total = n_velocities * nx * ny

    I_idx = Int[]
    J_idx = Int[]
    V_vals = Float64[]

    global_index(vel_idx, i, j) = ((j - 1) * nx + (i - 1)) * n_velocities + vel_idx

    for j in 1:ny
        for i in 1:nx
            for vel in 1:n_velocities
                row_idx = global_index(vel, i, j)
                ex, ey = e[vel]

                if ex == 0 && ey == 0
                    push!(I_idx, row_idx)
                    push!(J_idx, row_idx)
                    push!(V_vals, 0.0)
                    continue
                end

                if ex != 0
                    left_i = i == 1 ? nx : i - 1
                    right_i = i == nx ? 1 : i + 1

                    push!(I_idx, row_idx)
                    push!(J_idx, global_index(vel, left_i, j))
                    push!(V_vals, -ex / (2 * hx))

                    push!(I_idx, row_idx)
                    push!(J_idx, global_index(vel, right_i, j))
                    push!(V_vals, ex / (2 * hx))
                end

                if ey != 0
                    reflected_vel = opposite_velocity_index_D2Q9(vel)

                    if j > 1
                        push!(I_idx, row_idx)
                        push!(J_idx, global_index(vel, i, j - 1))
                        push!(V_vals, -ey / (2 * hy))
                    else
                        push!(I_idx, row_idx)
                        push!(J_idx, global_index(reflected_vel, i, j))
                        push!(V_vals, -ey / (2 * hy))
                    end

                    if j < ny
                        push!(I_idx, row_idx)
                        push!(J_idx, global_index(vel, i, j + 1))
                        push!(V_vals, ey / (2 * hy))
                    else
                        push!(I_idx, row_idx)
                        push!(J_idx, global_index(reflected_vel, i, j))
                        push!(V_vals, ey / (2 * hy))
                    end
                end
            end
        end
    end

    return sparse(I_idx, J_idx, V_vals, n_total, n_total), e
end

function tg_velocity_field(i, j, nx, ny, amplitude)
    x = 2 * π * (i - 1) / nx
    y = 2 * π * (j - 1) / ny
    ux = amplitude * sin(x) * cos(y)
    uy = -amplitude * cos(x) * sin(y)
    return ux, uy
end

function d2q9_equilibrium(rho, ux, uy, w_value, e_value, a_value, b_value, c_value, d_value)
    vcx = [velocity[1] for velocity in e_value]
    vcy = [velocity[2] for velocity in e_value]
    return cal_feq(rho, ux, uy, w_value, vcx, vcy, a_value, b_value, c_value, d_value)
end

function tg2d_initial_condition(nx, ny, amplitude, rho_value, w_value, e_value, a_value, b_value, c_value, d_value)
    phi_ini = zeros(9 * nx * ny)

    for j in 1:ny
        for i in 1:nx
            ux, uy = tg_velocity_field(i, j, nx, ny, amplitude)
            feq = d2q9_equilibrium(rho_value, ux, uy, w_value, e_value, a_value, b_value, c_value, d_value)
            start_idx = ((j - 1) * nx + (i - 1)) * 9 + 1
            phi_ini[start_idx:start_idx + 8] .= feq
        end
    end

    return phi_ini
end

function macroscopic_fields_from_state(phi, nx, ny, e_value)
    rho = zeros(nx, ny)
    ux = zeros(nx, ny)
    uy = zeros(nx, ny)

    for j in 1:ny
        for i in 1:nx
            start_idx = ((j - 1) * nx + (i - 1)) * 9 + 1
            fi = phi[start_idx:start_idx + 8]
            rho_local = sum(fi)
            rho[i, j] = rho_local
            ux[i, j] = sum(fi[k] * e_value[k][1] for k = 1:9) / rho_local
            uy[i, j] = sum(fi[k] * e_value[k][2] for k = 1:9) / rho_local
        end
    end

    return rho, ux, uy
end

function per_site_macroscopic(phi, Q, ngrid, e_val)
    n_time = size(phi, 2)
    rho_field = zeros(ngrid, n_time)
    umag_field = zeros(ngrid, n_time)

    for nt = 1:n_time
        f = reshape(phi[:, nt], Q, ngrid)
        for s = 1:ngrid
            rho_local = sum(f[:, s])
            rho_field[s, nt] = rho_local
            ux = sum(e_val[k][1] * f[k, s] for k = 1:Q) / max(rho_local, eps())
            uy = sum(e_val[k][2] * f[k, s] for k = 1:Q) / max(rho_local, eps())
            umag_field[s, nt] = sqrt(ux^2 + uy^2)
        end
    end

    return rho_field, umag_field
end

function domain_average_distribution_history(phiT, Q, ngrid)
    n_time = size(phiT, 2)
    avg_phi = zeros(Q, n_time)

    for nt = 1:n_time
        avg_phi[:, nt] = vec(mean(reshape(phiT[:, nt], Q, ngrid), dims=2))
    end

    return avg_phi
end

function density_error_history(phiT_ref, phiT_clbm, nx, ny, e_value)
    n_time = size(phiT_ref, 2)
    abs_err = zeros(n_time)

    for nt in 1:n_time
        rho_ref, _, _ = macroscopic_fields_from_state(phiT_ref[:, nt], nx, ny, e_value)
        rho_clbm, _, _ = macroscopic_fields_from_state(phiT_clbm[:, nt], nx, ny, e_value)
        abs_err[nt] = norm(vec(rho_clbm - rho_ref))
    end

    return abs_err
end

function velocity_error_history(phiT_ref, phiT_clbm, nx, ny, e_value)
    n_time = size(phiT_ref, 2)
    abs_err = zeros(n_time)
    rel_err = zeros(n_time)

    for nt in 1:n_time
        _, ux_ref, uy_ref = macroscopic_fields_from_state(phiT_ref[:, nt], nx, ny, e_value)
        _, ux_clbm, uy_clbm = macroscopic_fields_from_state(phiT_clbm[:, nt], nx, ny, e_value)

        err_vec = vcat(vec(ux_clbm - ux_ref), vec(uy_clbm - uy_ref))
        ref_vec = vcat(vec(ux_ref), vec(uy_ref))
        abs_err[nt] = norm(err_vec)
        rel_err[nt] = abs_err[nt] / max(norm(ref_vec), eps(Float64))
    end

    return abs_err, rel_err
end

function build_tg2d_diagnostics(phiT_ref, phiT_clbm, nx, ny, e_value)
    ngrid_local = nx * ny
    n_time_local = size(phiT_ref, 2)
    rho_ref_field, u_ref_field = per_site_macroscopic(phiT_ref, Q, ngrid_local, e_value)
    rho_clbm_field, u_clbm_field = per_site_macroscopic(phiT_clbm, Q, ngrid_local, e_value)

    rho_ref_mean = vec(mean(rho_ref_field, dims=1))
    rho_clbm_mean = vec(mean(rho_clbm_field, dims=1))
    u_ref_rms = vec([sqrt(mean(u_ref_field[:, nt].^2)) for nt = 1:n_time_local])
    u_clbm_rms = vec([sqrt(mean(u_clbm_field[:, nt].^2)) for nt = 1:n_time_local])

    avg_phi_ref = domain_average_distribution_history(phiT_ref, Q, ngrid_local)
    avg_phi_clbm = domain_average_distribution_history(phiT_clbm, Q, ngrid_local)
    avg_phi_abs_err = abs.(avg_phi_clbm .- avg_phi_ref)

    dist_abs_err = abs.(phiT_clbm .- phiT_ref)
    dist_rel_err = dist_abs_err ./ max.(abs.(phiT_ref), eps(Float64))
    profile_abs_max = vec(maximum(dist_abs_err, dims=1))
    profile_abs_l2 = vec([norm(dist_abs_err[:, nt]) for nt = 1:n_time_local])
    density_error_norm = density_error_history(phiT_ref, phiT_clbm, nx, ny, e_value)
    velocity_error_norm, velocity_rel_error_norm = velocity_error_history(phiT_ref, phiT_clbm, nx, ny, e_value)

    return (
        rho_ref_mean=rho_ref_mean,
        rho_clbm_mean=rho_clbm_mean,
        rho_mean_abs_err=abs.(rho_clbm_mean .- rho_ref_mean),
        u_ref_rms=u_ref_rms,
        u_clbm_rms=u_clbm_rms,
        u_rms_abs_err=abs.(u_clbm_rms .- u_ref_rms),
        dist_rel_err=dist_rel_err,
        density_error_norm=density_error_norm,
        velocity_error_norm=velocity_error_norm,
        velocity_rel_error_norm=velocity_rel_error_norm,
        profile_abs_max=profile_abs_max,
        profile_abs_l2=profile_abs_l2,
        avg_phi_ref=avg_phi_ref,
        avg_phi_clbm=avg_phi_clbm,
        avg_phi_abs_err=avg_phi_abs_err,
    )
end

function write_tg2d_snapshot_h5(filename, phi_ref, phi_clbm, nx, ny, e_value; time_step, time_value)
    rho_ref, ux_ref, uy_ref = macroscopic_fields_from_state(phi_ref, nx, ny, e_value)
    rho_clbm, ux_clbm, uy_clbm = macroscopic_fields_from_state(phi_clbm, nx, ny, e_value)
    umag_ref = sqrt.(ux_ref .^ 2 .+ uy_ref .^ 2)
    umag_clbm = sqrt.(ux_clbm .^ 2 .+ uy_clbm .^ 2)

    h5open(filename, "w") do file
        write(file, "time_step", time_step)
        write(file, "t", time_value)
        write(file, "phi_ref", reshape(phi_ref, Q, nx, ny))
        write(file, "phi_clbe", reshape(phi_clbm, Q, nx, ny))
        write(file, "abs_err_phi", reshape(abs.(phi_clbm .- phi_ref), Q, nx, ny))
        write(file, "rho_ref", rho_ref)
        write(file, "rho_clbe", rho_clbm)
        write(file, "abs_err_rho", abs.(rho_clbm .- rho_ref))
        write(file, "ux_ref", ux_ref)
        write(file, "uy_ref", uy_ref)
        write(file, "ux_clbe", ux_clbm)
        write(file, "uy_clbe", uy_clbm)
        write(file, "abs_err_ux", abs.(ux_clbm .- ux_ref))
        write(file, "abs_err_uy", abs.(uy_clbm .- uy_ref))
        write(file, "u_mag_ref", umag_ref)
        write(file, "u_mag_clbe", umag_clbm)
        write(file, "abs_err_u_mag", abs.(umag_clbm .- umag_ref))
    end
end

function save_tg2d_comparison_hdf5(result; output_dir="data/tg2d_clbe_comparison", snapshot_every=0)
    mkpath(output_dir)
    diagnostics = result.diagnostics
    snapshot_paths = String[]
    effective_snapshot_every = snapshot_every <= 0 ? max(1, result.local_n_time - 1) : snapshot_every
    snapshot_indices = sort(unique(vcat(1, collect(1:effective_snapshot_every:result.local_n_time), result.local_n_time)))

    for nt in snapshot_indices
        filename = joinpath(output_dir, "var" * @sprintf("%08d", nt - 1) * ".h5")
        write_tg2d_snapshot_h5(
            filename,
            result.phiT_ref[:, nt],
            result.phiT_clbm[:, nt],
            result.nx,
            result.ny,
            result.e_value;
            time_step=nt - 1,
            time_value=dt * (nt - 1),
        )
        push!(snapshot_paths, filename)
    end

    fn_ts = joinpath(output_dir, "time_series.h5")
    h5open(fn_ts, "w") do file
        write(file, "t", dt .* collect(0:result.local_n_time - 1))
        write(file, "density_ref_mean", diagnostics.rho_ref_mean)
        write(file, "density_clbe_mean", diagnostics.rho_clbm_mean)
        write(file, "density_mean_abs_err", diagnostics.rho_mean_abs_err)
        write(file, "velocity_ref_rms", diagnostics.u_ref_rms)
        write(file, "velocity_clbe_rms", diagnostics.u_clbm_rms)
        write(file, "velocity_rms_abs_err", diagnostics.u_rms_abs_err)
        write(file, "density_error_norm", diagnostics.density_error_norm)
        write(file, "velocity_error_norm", diagnostics.velocity_error_norm)
        write(file, "profile_abs_max", diagnostics.profile_abs_max)
        write(file, "profile_abs_l2", diagnostics.profile_abs_l2)
        write(file, "avg_phi_ref", diagnostics.avg_phi_ref)
        write(file, "avg_phi_clbe", diagnostics.avg_phi_clbm)
        write(file, "avg_phi_abs_err", diagnostics.avg_phi_abs_err)
        write(file, "case_label", result.case_label)
        write(file, "reference_model", String(result.reference_model))
        write(file, "nx", result.nx)
        write(file, "ny", result.ny)
        write(file, "truncation_order", result.local_truncation_order)
        write(file, "rho_value", result.rho_value)
        write(file, "amplitude", result.amplitude)
    end

    return (
        output_dir=output_dir,
        time_series_path=fn_ts,
        snapshot_paths=snapshot_paths,
    )
end

function plot_tg2d_comparison(phiT_ref, phiT_clbm, nx, ny, e_value, local_n_time, truncation_order; case_label="2D TG periodic test")
    _, ux_ref, uy_ref = macroscopic_fields_from_state(phiT_ref[:, end], nx, ny, e_value)
    _, ux_clbm, uy_clbm = macroscopic_fields_from_state(phiT_clbm[:, end], nx, ny, e_value)
    ux_err = ux_clbm - ux_ref
    uy_err = uy_clbm - uy_ref
    abs_err, _ = velocity_error_history(phiT_ref, phiT_clbm, nx, ny, e_value)

    close("all")
    figure(figsize=(14, 8))

    subplot(2, 3, 1)
    imshow(ux_ref', origin="lower", cmap="RdBu")
    colorbar()
    title(L"u_x^{\mathrm{LBM}}")

    subplot(2, 3, 2)
    imshow(ux_clbm', origin="lower", cmap="RdBu")
    colorbar()
    title(L"u_x^{\mathrm{CLBM}}")

    subplot(2, 3, 3)
    imshow(ux_err', origin="lower", cmap="RdBu")
    colorbar()
    title(L"u_x^{\mathrm{CLBM}} - u_x^{\mathrm{LBM}}")

    subplot(2, 3, 4)
    imshow(uy_ref', origin="lower", cmap="RdBu")
    colorbar()
    title(L"u_y^{\mathrm{LBM}}")

    subplot(2, 3, 5)
    imshow(uy_clbm', origin="lower", cmap="RdBu")
    colorbar()
    title(L"u_y^{\mathrm{CLBM}}")

    subplot(2, 3, 6)
    semilogy(1:local_n_time, abs_err, "-b", linewidth=1.8, label="absolute")
    xlabel("Time step")
    ylabel("Velocity error norm")
    title("Velocity error history")

    suptitle("$case_label, grid = $(nx)×$(ny), k = $truncation_order")
    tight_layout(rect=(0, 0, 1, 0.96))
    display(gcf())
    show()
end

function build_numerical_tg_reference(; nx, ny, amplitude, rho_value, tau_value, local_n_time, boundary_setup)
    case_label = boundary_setup ? "2D TG boundary-initialized test" : "2D TG periodic test"

    if boundary_setup
        println("Using the current D2Q9 LBM boundary-value TG initializer as the shared setup.")
        println("ℹ️  Reference LBM history uses the current top/bottom no-slip scheme from src/LBM/streaming.jl.")
        _, _, _, _, _, _, _, reference_phi_history = run_tg_d2q9_boundary_lbm(
            nx=nx,
            ny=ny,
            amplitude=amplitude,
            rho_value=rho_value,
            tau_value=tau_value,
            n_time=local_n_time,
            return_phi_history=true,
        )
    else
        println("Using the current pure numerical D2Q9 LBM as the periodic reference simulation.")
        _, _, _, _, _, reference_phi_history = run_tg_d2q9_lbm(
            nx=nx,
            ny=ny,
            amplitude=amplitude,
            rho_value=rho_value,
            tau_value=tau_value,
            n_time=local_n_time,
            l_noslipBC=false,
            return_phi_history=true,
        )
    end

    return (
        case_label=case_label,
        reference_initial_state=copy(reference_phi_history[:, 1]),
        reference_phi_history=reference_phi_history,
    )
end

function initialize_d2q9_tg_globals!(; nx, ny, rho_value, coeff_method, local_truncation_order)
    # QCFD convention: ngrid = LX * LY * LZ. D2Q9 TG flow is 2D, so LZ = 1.
    global LX = nx
    global LY = ny
    global LZ = 1
    global ngrid = LX * LY * LZ
    global Q = 9
    global D = 2
    global use_sparse = true
    global force_factor = 0.0
    global rho0 = rho_value
    global lTaylor = true
    global truncation_order = local_truncation_order
    global coeff_generation_method = coeff_method

    return ngrid
end

function build_direct_lbe_tg_reference(; nx, ny, amplitude, rho_value, local_n_time, boundary_setup, setup, direct_lbe_integrator=:euler)
    case_label = boundary_setup ? "2D TG direct n-point LBE boundary-aware test" : "2D TG direct n-point LBE periodic test"

    if boundary_setup
        println("Using direct n-point LBE with the same boundary-aware centered-difference D2Q9 streaming operator as CLBE.")
    else
        println("Using direct n-point LBE with the same periodic centered-difference D2Q9 streaming operator as CLBE.")
    end

    phi_ini = tg2d_initial_condition(
        nx,
        ny,
        amplitude,
        rho_value,
        setup.numeric_weights,
        setup.numeric_velocities,
        setup.a_val,
        setup.b_val,
        setup.c_val,
        setup.d_val,
    )

    hx = 1.0
    hy = 1.0
    S_lbm, _ = select_d2q9_streaming_operator(nx, ny, hx, hy; boundary_setup=boundary_setup)
    reference_phi_history = timeMarching_direct_LBE_ngrid(
        phi_ini,
        dt,
        local_n_time,
        setup.carleman_F1,
        setup.carleman_F2,
        setup.carleman_F3;
        S_lbm=S_lbm,
        integrator=direct_lbe_integrator,
    )

    return (
        case_label=case_label,
        reference_initial_state=phi_ini,
        reference_phi_history=reference_phi_history,
        S_lbm=S_lbm,
    )
end

function build_tg_reference(; nx, ny, amplitude, rho_value, tau_value, local_n_time, boundary_setup, reference_model, setup, direct_lbe_integrator=:euler)
    if reference_model == :direct_lbe
        return build_direct_lbe_tg_reference(
            nx=nx,
            ny=ny,
            amplitude=amplitude,
            rho_value=rho_value,
            local_n_time=local_n_time,
            boundary_setup=boundary_setup,
            setup=setup,
            direct_lbe_integrator=direct_lbe_integrator,
        )
    elseif reference_model == :lbm
        reference = build_numerical_tg_reference(
            nx=nx,
            ny=ny,
            amplitude=amplitude,
            rho_value=rho_value,
            tau_value=tau_value,
            local_n_time=local_n_time,
            boundary_setup=boundary_setup,
        )
        hx = 1.0
        hy = 1.0
        S_lbm, _ = select_d2q9_streaming_operator(nx, ny, hx, hy; boundary_setup=boundary_setup)
        return merge(reference, (S_lbm=S_lbm,))
    else
        error("Unsupported reference_model=$(reference_model). Supported values are :direct_lbe and :lbm.")
    end
end

function build_symbolic_carleman_setup(; rho_value, nspatial)
    println("Using symbolic D2Q9 LBM only to derive Carleman operators.")

    symbolic_weights, symbolic_velocities, numeric_weights, numeric_velocities, a_sym, b_sym, c_sym, d_sym, a_val, b_val, c_val, d_val = lbm_const_sym(Q_local=Q, D_local=D)
    symbolic_state, symbolic_collision, symbolic_velocity, symbolic_density = collision(Q, D, symbolic_weights, symbolic_velocities, rho_value, lTaylor, lorder2)
    carleman_F1, carleman_F2, carleman_F3 = get_coeff_LBM_Fi_ngrid(poly_order, Q, symbolic_state, symbolic_collision, tau_value, nspatial)

    return (
        symbolic_weights=symbolic_weights,
        symbolic_velocities=symbolic_velocities,
        numeric_weights=numeric_weights,
        numeric_velocities=numeric_velocities,
        a_sym=a_sym,
        b_sym=b_sym,
        c_sym=c_sym,
        d_sym=d_sym,
        a_val=a_val,
        b_val=b_val,
        c_val=c_val,
        d_val=d_val,
        symbolic_state=symbolic_state,
        symbolic_collision=symbolic_collision,
        symbolic_velocity=symbolic_velocity,
        symbolic_density=symbolic_density,
        carleman_F1=carleman_F1,
        carleman_F2=carleman_F2,
        carleman_F3=carleman_F3,
    )
end

function build_carleman_setup(; rho_value, nspatial, method=coeff_generation_method)
    numeric_weights, numeric_velocities = lbm_const_numerical(Q_local=Q, D_local=D)[1:2]

    if method == :symbolic
        symbolic_setup = build_symbolic_carleman_setup(rho_value=rho_value, nspatial=nspatial)
        return merge(symbolic_setup, (method=method,))
    elseif method == :numerical
        println("Using pure numerical coefficient generation for D2Q9 Carleman operators.")

        symbolic_weights, symbolic_velocities, _, _, a_sym, b_sym, c_sym, d_sym, a_val, b_val, c_val, d_val = lbm_const_sym(Q_local=Q, D_local=D)
        symbolic_state, symbolic_collision, symbolic_velocity, symbolic_density = collision(Q, D, symbolic_weights, symbolic_velocities, rho_value, lTaylor, lorder2)
        carleman_F1, carleman_F2, carleman_F3 = get_coeff_LBM_Fi_ngrid(
            poly_order,
            Q,
            symbolic_state,
            symbolic_collision,
            tau_value,
            nspatial;
            method=:numerical,
            w_value_input=numeric_weights,
            e_value_input=numeric_velocities,
            rho_value_input=rho_value,
            lTaylor_input=lTaylor,
            D_input=D,
        )

        return (
            method=method,
            symbolic_weights=symbolic_weights,
            symbolic_velocities=symbolic_velocities,
            numeric_weights=numeric_weights,
            numeric_velocities=numeric_velocities,
            a_sym=a_sym,
            b_sym=b_sym,
            c_sym=c_sym,
            d_sym=d_sym,
            a_val=a_val,
            b_val=b_val,
            c_val=c_val,
            d_val=d_val,
            symbolic_state=symbolic_state,
            symbolic_collision=symbolic_collision,
            symbolic_velocity=symbolic_velocity,
            symbolic_density=symbolic_density,
            carleman_F1=carleman_F1,
            carleman_F2=carleman_F2,
            carleman_F3=carleman_F3,
        )
    else
        error("Unknown coefficient-generation method: $method")
    end
end

function select_d2q9_streaming_operator(nx, ny, hx, hy; boundary_setup=false)
    if boundary_setup
        println("ℹ️  CLBM uses a boundary-aware centered-difference D2Q9 streaming operator with bounce-back-inspired wall coupling.")
        return streaming_operator_D2Q9_interleaved_boundary_aware(nx, ny, hx, hy)
    end

    return streaming_operator_D2Q9_interleaved_periodic(nx, ny, hx, hy)
end

function prepare_d2q9_carleman_runtime(; nx, ny, rho_value=1.0001, coeff_method=coeff_generation_method, local_truncation_order=truncation_order, boundary_setup=false, hx=1.0, hy=1.0)
    ngrid_local = initialize_d2q9_tg_globals!(
        nx=nx,
        ny=ny,
        rho_value=rho_value,
        coeff_method=coeff_method,
        local_truncation_order=local_truncation_order,
    )

    setup = build_carleman_setup(rho_value=rho_value, nspatial=ngrid_local, method=coeff_method)
    S_lbm, e_stream = select_d2q9_streaming_operator(nx, ny, hx, hy; boundary_setup=boundary_setup)

    global w_value = setup.numeric_weights
    global e_value = setup.numeric_velocities
    global F1_ngrid = setup.carleman_F1
    global F2_ngrid = setup.carleman_F2
    global F3_ngrid = setup.carleman_F3

    return (
        ngrid=ngrid_local,
        setup=setup,
        S_lbm=S_lbm,
        streaming_velocities=e_stream,
        w_value=w_value,
        e_value=e_value,
    )
end

function run_tg2d_clbe_comparison(; nx=3, ny=3, amplitude=0.05, rho_value=1.0001, local_n_time=n_time, boundary_setup=false, coeff_method=coeff_generation_method, local_truncation_order=truncation_order, reference_model=:direct_lbe, integrator=:euler, direct_lbe_integrator=:euler)
    if nx < 3 || ny < 3
        error("Use nx >= 3 and ny >= 3 for non-degenerate periodic centered-difference TG streaming.")
    end

    if local_truncation_order < poly_order
        error("local_truncation_order must satisfy local_truncation_order >= poly_order. Got local_truncation_order = $(local_truncation_order), poly_order = $(poly_order).")
    end

    runtime = prepare_d2q9_carleman_runtime(
        nx=nx,
        ny=ny,
        rho_value=rho_value,
        coeff_method=coeff_method,
        local_truncation_order=local_truncation_order,
        boundary_setup=boundary_setup,
    )
    ngrid_local = runtime.ngrid
    setup = runtime.setup

    reference = build_tg_reference(
        nx=nx,
        ny=ny,
        amplitude=amplitude,
        rho_value=rho_value,
        tau_value=tau_value,
        local_n_time=local_n_time,
        boundary_setup=boundary_setup,
        reference_model=reference_model,
        setup=setup,
        direct_lbe_integrator=direct_lbe_integrator,
    )

    phi_ini = reference.reference_initial_state
    phiT_ref = reference.reference_phi_history
    S_lbm = reference.S_lbm
    case_label = reference.case_label

    phiT_clbm, VT = timeMarching_state_CLBM_sparse(
        setup.symbolic_collision,
        setup.symbolic_state,
        tau_value,
        Q,
        truncation_order,
        dt,
        phi_ini,
        local_n_time;
        S_lbm=S_lbm,
        nspatial=ngrid_local,
        integrator=integrator,
    )

    dist_abs_err = abs.(phiT_clbm .- phiT_ref)
    dist_rel_err = dist_abs_err ./ max.(abs.(phiT_ref), eps(Float64))
    diagnostics = build_tg2d_diagnostics(phiT_ref, phiT_clbm, nx, ny, e_value)

    return (
        case_label=case_label,
        reference_model=reference_model,
        phi_ini=phi_ini,
        phiT_ref=phiT_ref,
        phiT_clbm=phiT_clbm,
        VT=VT,
        dist_abs_err=dist_abs_err,
        dist_rel_err=dist_rel_err,
        density_error_norm=diagnostics.density_error_norm,
        vel_abs_err=diagnostics.velocity_error_norm,
        vel_rel_err=diagnostics.velocity_rel_error_norm,
        diagnostics=diagnostics,
        e_value=e_value,
        setup=setup,
        S_lbm=S_lbm,
        nx=nx,
        ny=ny,
        local_n_time=local_n_time,
        local_truncation_order=local_truncation_order,
        rho_value=rho_value,
        amplitude=amplitude,
        saved_paths=nothing,
    )
end

function main(; nx=3, ny=3, amplitude=0.05, rho_value=1.0001, local_n_time=n_time, l_plot=false, boundary_setup=false, coeff_method=coeff_generation_method, local_truncation_order=truncation_order, reference_model=:direct_lbe, save_output=false, output_dir="data/tg2d_clbe_comparison", snapshot_every=0, integrator=:euler, direct_lbe_integrator=:euler)
    result = run_tg2d_clbe_comparison(
        nx=nx,
        ny=ny,
        amplitude=amplitude,
        rho_value=rho_value,
        local_n_time=local_n_time,
        boundary_setup=boundary_setup,
        coeff_method=coeff_method,
        local_truncation_order=local_truncation_order,
        reference_model=reference_model,
        integrator=integrator,
        direct_lbe_integrator=direct_lbe_integrator,
    )

    if save_output
        saved_paths = save_tg2d_comparison_hdf5(result; output_dir=output_dir, snapshot_every=snapshot_every)
        result = merge(result, (saved_paths=saved_paths,))
    end

    phiT_lbe = result.phiT_ref
    phiT_clbm = result.phiT_clbm
    VT = result.VT
    dist_abs_err = result.dist_abs_err
    dist_rel_err = result.dist_rel_err
    density_error_norm = result.density_error_norm
    vel_abs_err = result.vel_abs_err
    vel_rel_err = result.vel_rel_err
    case_label = result.case_label

    println("Running $(case_label) CLBE comparison")
    println("  grid = $(nx)×$(ny)")
    println("  Q = $Q, D = $D")
    println("  poly_order = $poly_order, truncation_order = $truncation_order")
    println("  n_time = $local_n_time")
    println("  reference_model = $(reference_model)")
    println("Max distribution absolute difference = ", maximum(dist_abs_err))
    println("Max density error norm = ", maximum(density_error_norm))
    println("Max velocity absolute error norm = ", maximum(vel_abs_err))
    if save_output
        println("Saved HDF5 outputs to ", result.saved_paths.output_dir)
        println("Saved time series to ", result.saved_paths.time_series_path)
    end

    if l_plot
        plot_tg2d_comparison(phiT_lbe, phiT_clbm, nx, ny, result.e_value, local_n_time, truncation_order; case_label=case_label)
    end

    return phiT_lbe, phiT_clbm, VT, dist_abs_err, dist_rel_err, vel_abs_err, vel_rel_err
end
