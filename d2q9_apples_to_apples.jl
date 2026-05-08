using Printf
using LinearAlgebra
using SparseArrays

ENV["QCFD_HOME"] = pwd()
ENV["QCFD_SRC"]  = pwd() * "/src/"

include("src/CLBE/clbe_tg2d_run.jl")

# D2Q9 apples-to-apples benchmark — mirrors the D1Q3 multigrid strategy:
#
#   * CLBM uses the assembled sparse (C - S) at truncation order k, marched by
#     explicit Euler with the repository's default dt = 0.1.
#   * Reference = direct semi-discrete n-point LBE, marched by explicit Euler
#     with the SAME dt, SAME F^(1), F^(2), F^(3), SAME centered-difference S.
#
# Both sides use rho_value = 1.0001 so the `rho_reference == 1.0 ?` shortcut in
# numerical_collision_rhs falls through to rho_local = sum(state), giving a
# well-posed polynomial decomposition (no frozen-ρ constant term).
#
# If the CLBM construction is correct, the only error this benchmark exposes
# is pure Carleman truncation at order k. The error should shrink with k and
# stay bounded over long nt, as it does in the D1Q3 multigrid test.

struct Case
    label::String
    nx::Int
    ny::Int
    k::Int
    amplitude::Float64
    n_time::Int
end

cases = Case[
    Case("nx3_k3_nt100", 3, 3, 3, 0.02, 100),
    Case("nx3_k4_nt100", 3, 3, 4, 0.02, 100),
    Case("nx4_k3_nt100", 4, 4, 3, 0.02, 100),
    Case("nx6_k3_nt100", 6, 6, 3, 0.02, 100),
]

function run_case(case::Case)
    t0 = time()

    result = run_tg2d_clbe_comparison(
        nx=case.nx,
        ny=case.ny,
        amplitude=case.amplitude,
        rho_value=1.0001,
        local_n_time=case.n_time,
        boundary_setup=false,
        coeff_method=:numerical,
        local_truncation_order=case.k,
        reference_model=:direct_lbe,
    )
    save_dir = joinpath("data", "d2q9_apples_to_apples", case.label)
    saved_paths = save_tg2d_comparison_hdf5(result; output_dir=save_dir, snapshot_every=max(1, case.n_time ÷ 10))

    elapsed = time() - t0
    abs_err = result.dist_abs_err
    avg_abs_err = result.diagnostics.avg_phi_abs_err

    return (
        case = case,
        elapsed = elapsed,
        abs_err = abs_err,
        avg_abs_err = avg_abs_err,
        density_error_norm = result.density_error_norm,
        velocity_error_norm = result.vel_abs_err,
        saved_paths = saved_paths,
    )
end

using Statistics

summaries = String[]

for case in cases
    println("\n" * "="^72)
    @printf("Running %s: nx=ny=%d, k=%d, A=%.3f, n_time=%d\n",
            case.label, case.nx, case.k, case.amplitude, case.n_time)
    println("="^72)
    GC.gc()
    r = run_case(case)

    println()
    @printf "Wall time: %.1f s\n" r.elapsed
    @printf "Per-site distribution error (|φ_CLBM − φ_direct|):\n"
    @printf "  overall max |Δφ|          = %.3e\n" maximum(r.abs_err)
    @printf "  final-step max |Δφ|       = %.3e\n" maximum(r.abs_err[:, end])
    @printf "Domain-averaged error (|<φ_CLBM> − <φ_direct>|):\n"
    @printf "  overall max |Δ<φ>|        = %.3e\n" maximum(r.avg_abs_err)
    @printf "  final-step max |Δ<φ>|     = %.3e\n" maximum(r.avg_abs_err[:, end])
    @printf "  max density error norm    = %.3e\n" maximum(r.density_error_norm)
    @printf "  max velocity error norm   = %.3e\n" maximum(r.velocity_error_norm)
    @printf "  per-component (max abs / final abs):\n"
    for m = 1:9
        @printf "    f_%d: max = %.3e   final = %.3e\n" m maximum(r.avg_abs_err[m, :]) r.avg_abs_err[m, end]
    end
    println("  saved time series: $(r.saved_paths.time_series_path)")

    push!(summaries, @sprintf("%-16s  nx=%2d  k=%d  nt=%3d  |Δφ|_max=%.3e  |Δφ|_final=%.3e  |Δ<φ>|_max=%.3e",
        case.label, case.nx, case.k, case.n_time,
        maximum(r.abs_err), maximum(r.abs_err[:, end]), maximum(r.avg_abs_err)))
end

println("\n\n" * "="^88)
println("D2Q9 CLBM vs direct n-point LBE  (apples-to-apples, same F/S/dt, rho0=1.0001)")
println("="^88)
for s in summaries; println(s); end

# Save summary
open("data/d2q9_apples_to_apples_summary.txt", "w") do io
    println(io, "D2Q9 CLBM vs direct n-point LBE benchmark")
    println(io, "rho0=1.0001, lTaylor=true, coeff_method=:numerical, dt=0.1, periodic TG")
    println(io, "Self-consistency: both sides use the same F^(1/2/3), same S, same dt.")
    println(io, "This mirrors the D1Q3 plot_multigrid_domain_average benchmark.")
    println(io)
    for s in summaries; println(io, s); end
end
println("\nSaved summary to data/d2q9_apples_to_apples_summary.txt")
