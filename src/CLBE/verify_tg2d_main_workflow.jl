using Test
using Statistics
using Printf

# Direct end-to-end benchmark of the top-level D2Q9 TG driver workflow.
#
# This script benchmarks src/CLBE/clbe_tg2d_run.jl::main(...) itself, rather than
# reconstructing operator pieces manually. The goal is to exercise:
#   - numerical reference generation,
#   - Carleman setup,
#   - streaming-operator selection,
#   - sparse CLBM time marching,
#   - returned error arrays and consistency of reported diagnostics.
#
# The cases are intentionally small so they remain feasible for repeated checks.

include("clbe_tg2d_run.jl")

struct TG2DCase
    label::String
    nx::Int
    ny::Int
    amplitude::Float64
    rho_value::Float64
    local_n_time::Int
    boundary_setup::Bool
    coeff_method::Symbol
    local_truncation_order::Int
end

function summarize_case(case::TG2DCase)
    phiT_lbe, phiT_clbm, VT, dist_abs_err, dist_rel_err, vel_abs_err, vel_rel_err = main(
        nx=case.nx,
        ny=case.ny,
        amplitude=case.amplitude,
        rho_value=case.rho_value,
        local_n_time=case.local_n_time,
        l_plot=false,
        boundary_setup=case.boundary_setup,
        coeff_method=case.coeff_method,
        local_truncation_order=case.local_truncation_order,
    )

    recomputed_dist_abs = abs.(phiT_clbm .- phiT_lbe)
    recomputed_dist_rel = recomputed_dist_abs ./ max.(abs.(phiT_lbe), eps(Float64))
    recomputed_vel_abs, recomputed_vel_rel = velocity_error_history(phiT_lbe, phiT_clbm, case.nx, case.ny, e_value)

    first_step_dist_abs = size(dist_abs_err, 2) >= 2 ? maximum(dist_abs_err[:, 2]) : 0.0
    first_step_dist_rel = size(dist_rel_err, 2) >= 2 ? maximum(dist_rel_err[:, 2]) : 0.0
    first_step_vel_abs = length(vel_abs_err) >= 2 ? vel_abs_err[2] : 0.0
    first_step_vel_rel = length(vel_rel_err) >= 2 ? vel_rel_err[2] : 0.0

    final_step_dist_abs = maximum(dist_abs_err[:, end])
    final_step_dist_rel = maximum(dist_rel_err[:, end])
    final_step_vel_abs = vel_abs_err[end]
    final_step_vel_rel = vel_rel_err[end]

    return (
        phiT_lbe=phiT_lbe,
        phiT_clbm=phiT_clbm,
        VT=VT,
        dist_abs_err=dist_abs_err,
        dist_rel_err=dist_rel_err,
        vel_abs_err=vel_abs_err,
        vel_rel_err=vel_rel_err,
        recomputed_dist_abs=recomputed_dist_abs,
        recomputed_dist_rel=recomputed_dist_rel,
        recomputed_vel_abs=recomputed_vel_abs,
        recomputed_vel_rel=recomputed_vel_rel,
        max_dist_abs=maximum(dist_abs_err),
        max_dist_rel=maximum(dist_rel_err),
        max_vel_abs=maximum(vel_abs_err),
        max_vel_rel=maximum(vel_rel_err),
        first_step_dist_abs=first_step_dist_abs,
        first_step_dist_rel=first_step_dist_rel,
        first_step_vel_abs=first_step_vel_abs,
        first_step_vel_rel=first_step_vel_rel,
        final_step_dist_abs=final_step_dist_abs,
        final_step_dist_rel=final_step_dist_rel,
        final_step_vel_abs=final_step_vel_abs,
        final_step_vel_rel=final_step_vel_rel,
    )
end

function render_summary(case::TG2DCase, metrics)
    mode = case.boundary_setup ? "boundary-aware" : "periodic"
    return join([
        "Case: $(case.label)",
        "  mode = $mode",
        "  grid = $(case.nx)x$(case.ny)",
        "  truncation_order = $(case.local_truncation_order)",
        "  n_time = $(case.local_n_time)",
        "  amplitude = $(case.amplitude)",
        @sprintf("  first_step max dist abs err = %.6e", metrics.first_step_dist_abs),
        @sprintf("  first_step max dist rel err = %.6e", metrics.first_step_dist_rel),
        @sprintf("  first_step vel abs err     = %.6e", metrics.first_step_vel_abs),
        @sprintf("  first_step vel rel err     = %.6e", metrics.first_step_vel_rel),
        @sprintf("  final_step max dist abs err = %.6e", metrics.final_step_dist_abs),
        @sprintf("  final_step max dist rel err = %.6e", metrics.final_step_dist_rel),
        @sprintf("  final_step vel abs err      = %.6e", metrics.final_step_vel_abs),
        @sprintf("  final_step vel rel err      = %.6e", metrics.final_step_vel_rel),
        @sprintf("  overall max dist abs err = %.6e", metrics.max_dist_abs),
        @sprintf("  overall max dist rel err = %.6e", metrics.max_dist_rel),
        @sprintf("  overall max vel abs err  = %.6e", metrics.max_vel_abs),
        @sprintf("  overall max vel rel err  = %.6e", metrics.max_vel_rel),
    ], "\n")
end

cases = [
    TG2DCase("periodic_smoke", 3, 3, 0.05, 1.0, 5, false, coeff_generation_method, 3),
    TG2DCase("periodic_low_amplitude", 3, 3, 0.02, 1.0, 5, false, coeff_generation_method, 3),
    TG2DCase("boundary_smoke", 3, 3, 0.05, 1.0, 5, true, coeff_generation_method, 3),
]

summaries = String[]
periodic_metrics = Dict{String, Any}()

@testset "Direct TG2D main workflow benchmark" begin
    for case in cases
        @testset "$(case.label)" begin
            metrics = summarize_case(case)
            push!(summaries, render_summary(case, metrics))

            if !case.boundary_setup
                periodic_metrics[case.label] = metrics
            end

            @test size(metrics.phiT_lbe) == size(metrics.phiT_clbm)
            @test size(metrics.phiT_lbe, 2) == case.local_n_time
            @test size(metrics.VT, 2) == case.local_n_time
            @test size(metrics.VT, 1) >= size(metrics.phiT_clbm, 1)

            @test metrics.VT[1:size(metrics.phiT_clbm, 1), :] ≈ metrics.phiT_clbm atol=1e-12 rtol=1e-12

            @test metrics.phiT_lbe[:, 1] ≈ metrics.phiT_clbm[:, 1] atol=1e-12 rtol=1e-12
            @test metrics.dist_abs_err[:, 1] ≈ zeros(size(metrics.dist_abs_err, 1)) atol=1e-14 rtol=1e-14
            @test metrics.vel_abs_err[1] ≈ 0.0 atol=1e-14 rtol=1e-14
            @test metrics.vel_rel_err[1] ≈ 0.0 atol=1e-14 rtol=1e-14

            @test metrics.dist_abs_err ≈ metrics.recomputed_dist_abs atol=1e-12 rtol=1e-12
            @test metrics.dist_rel_err ≈ metrics.recomputed_dist_rel atol=1e-12 rtol=1e-12
            @test metrics.vel_abs_err ≈ metrics.recomputed_vel_abs atol=1e-12 rtol=1e-12
            @test metrics.vel_rel_err ≈ metrics.recomputed_vel_rel atol=1e-12 rtol=1e-12

            @test all(isfinite.(metrics.phiT_lbe))
            @test all(isfinite.(metrics.phiT_clbm))
            @test all(isfinite.(metrics.VT))
            @test all(isfinite.(metrics.dist_abs_err))
            @test all(isfinite.(metrics.dist_rel_err))
            @test all(isfinite.(metrics.vel_abs_err))
            @test all(isfinite.(metrics.vel_rel_err))

            # Hard workflow checks only: keep these structural and consistency-based.
            # The numerical error magnitudes are recorded in the summary file and should
            # be interpreted separately from pass/fail workflow integrity.
            @test metrics.max_dist_abs >= 0.0
            @test metrics.max_vel_abs >= 0.0
            @test metrics.max_vel_rel >= 0.0
        end
    end

    @testset "Periodic case cross-checks" begin
        smoke = periodic_metrics["periodic_smoke"]
        low_amp = periodic_metrics["periodic_low_amplitude"]

        @test smoke.first_step_dist_abs >= 0.0
        @test low_amp.first_step_dist_abs >= 0.0
        @test smoke.max_dist_abs >= low_amp.max_dist_abs
        @test smoke.max_vel_abs >= low_amp.max_vel_abs
    end
end

summary_path = joinpath("data", "verify_tg2d_main_workflow_summary.txt")
mkpath(dirname(summary_path))
open(summary_path, "w") do io
    println(io, "Direct TG2D main workflow benchmark summary")
    println(io, "==========================================")
    println(io)
    for block in summaries
        println(io, block)
        println(io)
    end
end

println("\nSaved TG2D main-workflow benchmark summary to: $summary_path")
for block in summaries
    println(block)
    println()
end

println("🎉 Direct TG2D main workflow benchmark completed successfully!")