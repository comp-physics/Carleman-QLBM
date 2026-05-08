using Printf
using SparseArrays
using LinearAlgebra

ENV["QCFD_HOME"] = pwd()
ENV["QCFD_SRC"]  = pwd() * "/src/"

include("src/CLBE/clbe_tg2d_run.jl")

# D2Q9 convergence sweep AFTER:
#   1) reverting the rho_local change (original logic)
#   2) changing lTaylor=false -> true in clbe_tg2d_run.jl::main() (the Carleman-
#      valid setting, as you noted).
# Uses dt=1 to match LBM lattice-time units, so CLBM step n and LBM iteration n
# are at the same physical time.

struct Case
    label::String
    nx::Int
    ny::Int
    k::Int
    amplitude::Float64
    n_time::Int
end

cases = Case[
    Case("nx3_k3_dt1",  3, 3, 3, 0.02, 10),
    Case("nx3_k4_dt1",  3, 3, 4, 0.02, 10),
    Case("nx4_k3_dt1",  4, 4, 3, 0.02, 10),
    Case("nx6_k3_dt1",  6, 6, 3, 0.02, 10),
]

results = Dict{String, NamedTuple}()

for case in cases
    println("\n" * "="^72)
    @printf("Running %s: nx=ny=%d, k=%d, A=%.3f, n_time=%d, dt=1\n",
            case.label, case.nx, case.k, case.amplitude, case.n_time)
    println("="^72)

    GC.gc()
    global dt = 1.0
    t0 = time()
    local ok, phiT_lbe, phiT_clbm, VT, dist_abs_err, dist_rel_err, vel_abs_err, vel_rel_err
    ok = false
    try
        phiT_lbe, phiT_clbm, VT, dist_abs_err, dist_rel_err, vel_abs_err, vel_rel_err = main(
            nx=case.nx, ny=case.ny, amplitude=case.amplitude, rho_value=1.0,
            local_n_time=case.n_time, l_plot=false, boundary_setup=false,
            coeff_method=:numerical, local_truncation_order=case.k,
        )
        ok = true
    catch err
        @warn "Case $(case.label) failed" exception=(err, catch_backtrace())
    end

    elapsed = time() - t0
    if ok
        # Per-step table (since CLBM time = LBM iteration at dt=1)
        nt_final = size(phiT_clbm, 2)
        println("\nPer-step comparison at dt=1 (matching physical time):")
        @printf "%4s  %12s  %12s  %12s  %12s\n" "step" "|φ_LBM|∞" "|φ_CLBM|∞" "|Δφ|∞" "|Δφ|_rel"
        for nt in [1, 2, 3, 5, min(10, nt_final)]
            if nt > nt_final; continue; end
            lbm_max = maximum(abs.(phiT_lbe[:, nt]))
            clbm_max = maximum(abs.(phiT_clbm[:, nt]))
            dphi = maximum(dist_abs_err[:, nt])
            @printf "%4d  %12.3e  %12.3e  %12.3e  %12.3e\n" (nt-1) lbm_max clbm_max dphi (dphi/max(lbm_max,eps()))
        end

        results[case.label] = (
            case = case,
            wall_s = elapsed,
            first_step_dist = maximum(dist_abs_err[:, 2]),
            first_step_vel  = vel_abs_err[2],
            overall_dist    = maximum(dist_abs_err),
            overall_vel     = maximum(vel_abs_err),
            final_dist      = maximum(dist_abs_err[:, end]),
            final_vel       = vel_abs_err[end],
            ok = true,
        )
    else
        results[case.label] = (case=case, wall_s=elapsed, ok=false)
    end
end

println("\n" * "="^88)
println("SUMMARY (D2Q9 periodic TG, dt=1, lTaylor=true, numerical, apples-to-apples LBM compare)")
println("="^88)
@printf "%-14s %3s %3s %6s %8s %12s %12s %12s %12s\n" "case" "nx" "k" "ntime" "wall[s]" "fs|Δf|" "fs|Δu|" "ov|Δf|" "ov|Δu|"
for case in cases
    r = results[case.label]
    if get(r, :ok, false)
        @printf "%-14s %3d %3d %6d %8.1f %12.3e %12.3e %12.3e %12.3e\n" case.label case.nx case.k case.n_time r.wall_s r.first_step_dist r.first_step_vel r.overall_dist r.overall_vel
    else
        @printf "%-14s %3d %3d %6d %8.1f  FAILED\n" case.label case.nx case.k case.n_time r.wall_s
    end
end
