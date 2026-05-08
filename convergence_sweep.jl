using Printf
using SparseArrays
using LinearAlgebra

ENV["QCFD_HOME"] = pwd()
ENV["QCFD_SRC"]  = pwd() * "/src/"

include("src/CLBE/clbe_tg2d_run.jl")

# ---------------------------------------------------------------------------
# Convergence sweep for D2Q9 CLBE vs reference LBM on the periodic TG vortex.
# For every case, the reference LBM is run on the SAME nx, ny grid inside
# clbe_tg2d_run.jl::main(...) → build_numerical_tg_reference(...), so the
# comparison is apples-to-apples by construction.
# ---------------------------------------------------------------------------

struct Case
    label::String
    nx::Int
    ny::Int
    k::Int
    amplitude::Float64
    n_time::Int
end

cases = Case[
    Case("nx3_k3_baseline",   3, 3, 3, 0.02, 10),
    Case("nx3_k4_truncation", 3, 3, 4, 0.02, 10),
    Case("nx4_k3",            4, 4, 3, 0.02, 10),
    Case("nx6_k3",            6, 6, 3, 0.02, 10),
    Case("nx8_k3",            8, 8, 3, 0.02, 10),
]

results = Dict{String, NamedTuple}()

for case in cases
    println("\n" * "="^72)
    @printf("Running case %s: nx=ny=%d, k=%d, A=%.3f, n_time=%d\n",
            case.label, case.nx, case.k, case.amplitude, case.n_time)
    println("="^72)

    GC.gc()
    t_start = time()

    local ok, phiT_lbe, phiT_clbm, VT, dist_abs_err, dist_rel_err, vel_abs_err, vel_rel_err
    ok = false
    try
        phiT_lbe, phiT_clbm, VT, dist_abs_err, dist_rel_err, vel_abs_err, vel_rel_err = main(
            nx               = case.nx,
            ny               = case.ny,
            amplitude        = case.amplitude,
            rho_value        = 1.0,
            local_n_time     = case.n_time,
            l_plot           = false,
            boundary_setup   = false,
            coeff_method     = :numerical,
            local_truncation_order = case.k,
        )
        ok = true
    catch err
        @warn "Case $(case.label) failed" exception=(err, catch_backtrace())
    end

    t_end = time()

    if ok
        first_step_dist_abs = size(dist_abs_err, 2) >= 2 ? maximum(dist_abs_err[:, 2]) : NaN
        first_step_vel_abs  = length(vel_abs_err) >= 2 ? vel_abs_err[2] : NaN
        final_dist_abs      = maximum(dist_abs_err[:, end])
        final_vel_abs       = vel_abs_err[end]
        overall_dist_abs    = maximum(dist_abs_err)
        overall_vel_abs     = maximum(vel_abs_err)
        ref_vel_final = begin
            _, ux_ref, uy_ref = macroscopic_fields_from_state(phiT_lbe[:, end], case.nx, case.ny, e_value)
            sqrt(sum(ux_ref .^ 2) + sum(uy_ref .^ 2))
        end
        results[case.label] = (
            case = case,
            wall_time_s = t_end - t_start,
            first_step_dist_abs = first_step_dist_abs,
            first_step_vel_abs  = first_step_vel_abs,
            final_dist_abs      = final_dist_abs,
            final_vel_abs       = final_vel_abs,
            overall_dist_abs    = overall_dist_abs,
            overall_vel_abs     = overall_vel_abs,
            ref_vel_norm_final  = ref_vel_final,
            succeeded = true,
        )
    else
        results[case.label] = (
            case = case,
            wall_time_s = t_end - t_start,
            first_step_dist_abs = NaN,
            first_step_vel_abs  = NaN,
            final_dist_abs      = NaN,
            final_vel_abs       = NaN,
            overall_dist_abs    = NaN,
            overall_vel_abs     = NaN,
            ref_vel_norm_final  = NaN,
            succeeded = false,
        )
    end
end

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

println("\n\n" * "="^80)
println("CONVERGENCE SWEEP SUMMARY (D2Q9 periodic TG, numerical coefficient path)")
println("="^80)
@printf("%-20s %4s %3s %6s %7s %10s %10s %10s %10s %10s\n",
    "case", "nx", "k", "ntime", "wall[s]",
    "fs|Δf|", "fs|Δu|", "ov|Δf|", "ov|Δu|", "‖u_ref‖")
for case in cases
    r = results[case.label]
    if r.succeeded
        @printf("%-20s %4d %3d %6d %7.1f %10.3e %10.3e %10.3e %10.3e %10.3e\n",
            case.label, case.nx, case.k, case.n_time, r.wall_time_s,
            r.first_step_dist_abs, r.first_step_vel_abs,
            r.overall_dist_abs,    r.overall_vel_abs,
            r.ref_vel_norm_final)
    else
        @printf("%-20s %4d %3d %6d %7.1f %s\n",
            case.label, case.nx, case.k, case.n_time, r.wall_time_s,
            "FAILED")
    end
end

# Persist summary
summary_path = "data/convergence_sweep_summary.txt"
open(summary_path, "w") do io
    println(io, "Convergence sweep summary")
    println(io, "=========================")
    println(io, "D2Q9 periodic Taylor-Green vortex. CLBM vs pure LBM reference on")
    println(io, "the same nx, ny grid (apples-to-apples). coeff_method=:numerical.")
    println(io, "Columns:")
    println(io, "  fs|Δf|  = max distribution abs err at first time step")
    println(io, "  fs|Δu|  = velocity abs err norm at first time step")
    println(io, "  ov|Δf|  = overall max distribution abs err over history")
    println(io, "  ov|Δu|  = overall max velocity abs err norm over history")
    println(io, "  ‖u_ref‖ = L2 norm of reference LBM velocity at final step")
    println(io)
    @printf(io, "%-20s %4s %3s %6s %7s %10s %10s %10s %10s %10s\n",
        "case", "nx", "k", "ntime", "wall[s]",
        "fs|Δf|", "fs|Δu|", "ov|Δf|", "ov|Δu|", "‖u_ref‖")
    for case in cases
        r = results[case.label]
        if r.succeeded
            @printf(io, "%-20s %4d %3d %6d %7.1f %10.3e %10.3e %10.3e %10.3e %10.3e\n",
                case.label, case.nx, case.k, case.n_time, r.wall_time_s,
                r.first_step_dist_abs, r.first_step_vel_abs,
                r.overall_dist_abs, r.overall_vel_abs,
                r.ref_vel_norm_final)
        else
            @printf(io, "%-20s %4d %3d %6d %7.1f %s\n",
                case.label, case.nx, case.k, case.n_time, r.wall_time_s,
                "FAILED")
        end
    end
end

println("\nSaved convergence sweep summary to: ", summary_path)
