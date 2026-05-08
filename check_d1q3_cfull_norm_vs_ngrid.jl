using Printf
using LinearAlgebra
using SparseArrays
using Random
using Dates
using HDF5
using PyPlot

if !haskey(ENV, "QCFD_HOME")
    ENV["QCFD_HOME"] = pwd()
end
if !haskey(ENV, "QCFD_SRC")
    ENV["QCFD_SRC"] = joinpath(pwd(), "src") * "/"
end

include("src/CLBE/clbe_multigrid_run.jl")
include("src/CLBE/cfull_norm_sweep_utils.jl")

const CFSU = CFullNormSweepUtils

const NGRID_SWEEP = [1, 3, 6, 12, 24, 48]
const DEFAULT_SWEEP_MODE = :explicit
const DEFAULT_START_NGRID = 3
const DEFAULT_STEP_NGRID = 3
const DEFAULT_END_NGRID = 12
const DEFAULT_INCREASE_FACTOR = 2
const DEFAULT_NUM_NGRID = 4
const MAX_SPECTRAL_DIM = 3_500_000

function tex_output_dir()
    return dirname(CFSU.figure_output_path("d1q3_cfull_norm_vs_ngrid.pdf"))
end

function summary_output_path()
    return CFSU.data_output_path("d1q3_cfull_norm_vs_ngrid_summary.txt")
end

function h5_output_path()
    return CFSU.data_output_path("d1q3_cfull_norm_vs_ngrid.h5")
end

function plot_output_path()
    return CFSU.figure_output_path("d1q3_cfull_norm_vs_ngrid.pdf")
end

function lifted_dimension(Q, truncation_order, ngrid)
    return sum((Q * ngrid)^k for k in 1:truncation_order)
end

function power_iteration_spectral_norm(A::SparseMatrixCSC; maxiter::Int=80, tol::Float64=1e-8, seed::Int=1234)
    return CFSU.power_iteration_spectral_norm(A; maxiter=maxiter, tol=tol, seed=seed)
end

function generate_ngrid_sweep(; start_ngrid::Int=DEFAULT_START_NGRID,
    step_ngrid::Int=DEFAULT_STEP_NGRID,
    end_ngrid::Int=DEFAULT_END_NGRID,
    increase_factor::Int=DEFAULT_INCREASE_FACTOR,
    num_ngrid::Int=DEFAULT_NUM_NGRID,
    sweep_mode::Symbol=:arithmetic)

    if sweep_mode == :arithmetic
        if step_ngrid <= 0
            error("step_ngrid must be positive, got $step_ngrid")
        end
        if end_ngrid < start_ngrid
            error("end_ngrid must be >= start_ngrid, got start_ngrid=$start_ngrid and end_ngrid=$end_ngrid")
        end
        return collect(start_ngrid:step_ngrid:end_ngrid)
    elseif sweep_mode == :geometric
        if increase_factor < 2
            error("increase_factor must be >= 2 for geometric sweeps, got $increase_factor")
        end
        if num_ngrid <= 0
            error("num_ngrid must be positive, got $num_ngrid")
        end
        return [start_ngrid * increase_factor^(i - 1) for i in 1:num_ngrid]
    else
        error("Unsupported generated sweep_mode=$sweep_mode. Use :arithmetic or :geometric")
    end
end

function resolve_ngrid_sweep(; ngrid_sweep=NGRID_SWEEP,
    start_ngrid::Int=DEFAULT_START_NGRID,
    step_ngrid::Int=DEFAULT_STEP_NGRID,
    end_ngrid::Int=DEFAULT_END_NGRID,
    increase_factor::Int=DEFAULT_INCREASE_FACTOR,
    num_ngrid::Int=DEFAULT_NUM_NGRID,
    sweep_mode::Symbol=DEFAULT_SWEEP_MODE)

    if sweep_mode == :explicit
        return collect(ngrid_sweep)
    elseif sweep_mode == :arithmetic || sweep_mode == :geometric
        return generate_ngrid_sweep(
            start_ngrid=start_ngrid,
            step_ngrid=step_ngrid,
            end_ngrid=end_ngrid,
            increase_factor=increase_factor,
            num_ngrid=num_ngrid,
            sweep_mode=sweep_mode,
        )
    else
        error("Unsupported sweep_mode=$sweep_mode. Use :explicit, :arithmetic, or :geometric")
    end
end

function build_full_generator(ngrid_val; truncation_order_val=3, coeff_method=:numerical)
    core_cfg, case_cfg = build_d1q3_multigrid_configs(
        comparison_ngrid=ngrid_val,
        local_use_sparse=true,
        local_n_time=1,
        local_truncation_order=truncation_order_val,
        coeff_method=coeff_method,
        initial_condition=:legacy,
        dt_override=1.0,
    )

    runtime = prepare_d1q3_runtime(core_cfg, case_cfg)
    S_lbm = ngrid_val > 1 ? streaming_operator_D1Q3_interleaved(ngrid_val, 1.0)[1] : nothing

    C_full, bt, _ = build_full_clbe_generator_sparse(
        runtime.omega,
        runtime.f,
        core_cfg.tau_value,
        core_cfg.Q,
        case_cfg.truncation_order,
        core_cfg.poly_order,
        core_cfg.force_factor,
        runtime.w_value,
        runtime.e_value;
        S_lbm=S_lbm,
        nspatial=ngrid_val,
    )

    dropzeros!(C_full)
    return C_full, bt, core_cfg, case_cfg
end

function write_summary(path, results)
    headers = ["ngrid", "dim", "nnz", "spectral_est", "inf_norm", "conv", "iters", "status"]
    rows = [[r.ngrid, r.dim, r.nnz, r.spectral_est, r.inf_norm, string(r.converged), r.iterations, r.status] for r in results]
    CFSU.write_summary_table(path;
        title="D1Q3 full CLBE generator norm sweep",
        intro_lines=[
            "Target matrix: C_full = C - S",
            "Norm columns: inf_norm is exact; spectral_est is a power-iteration estimate of the spectral/operator 2-norm when computed",
        ],
        headers=headers,
        rows=rows,
        left_align=[false, false, false, false, false, true, false, true],
    )
end

function save_results_h5(path, results)
    metadata = Dict(
        "target_matrix" => "C_full = C - S",
        "max_spectral_dim" => Int(MAX_SPECTRAL_DIM),
        "Q" => isempty(results) ? 3 : Int(results[1].Q),
        "truncation_order" => isempty(results) ? 3 : Int(results[1].truncation_order),
        "sweep_kind" => "ngrid",
    )
    CFSU.write_results_h5(path, results; metadata=metadata)
end

function load_results_h5(path)
    return CFSU.load_results_h5(path)
end

function plot_results(results, output_pdf)
    CFSU.plot_norm_results(results, output_pdf;
        x_values=[r.ngrid for r in results],
        xlabel_text="ngrid",
        title_text="D1Q3 full CLBE generator norm vs ngrid",
        include_inf=true,
        include_spectral=true,
        xscale_base=2,
    )
end

function plot_results_from_h5(; h5_path=h5_output_path(), output_pdf=plot_output_path())
    results = load_results_h5(h5_path)
    plot_results(results, output_pdf)
    return results
end

function main(; ngrid_sweep=NGRID_SWEEP,
    start_ngrid::Int=DEFAULT_START_NGRID,
    step_ngrid::Int=DEFAULT_STEP_NGRID,
    end_ngrid::Int=DEFAULT_END_NGRID,
    increase_factor::Int=DEFAULT_INCREASE_FACTOR,
    num_ngrid::Int=DEFAULT_NUM_NGRID,
    sweep_mode::Symbol=DEFAULT_SWEEP_MODE,
    plot::Bool=(get(ENV, "D1Q3_CFULL_MAKE_PLOT", "1") == "1"),
    load_h5_only::Bool=(get(ENV, "D1Q3_CFULL_LOAD_H5_ONLY", "0") == "1"))

    resolved_ngrid_sweep = resolve_ngrid_sweep(
        ngrid_sweep=ngrid_sweep,
        start_ngrid=start_ngrid,
        step_ngrid=step_ngrid,
        end_ngrid=end_ngrid,
        increase_factor=increase_factor,
        num_ngrid=num_ngrid,
        sweep_mode=sweep_mode,
    )

    output_pdf = plot_output_path()
    summary_path = summary_output_path()
    h5_path = h5_output_path()

    if load_h5_only
        if !isfile(h5_path)
            error("Cannot load cached D1Q3 norm data because HDF5 file does not exist: $h5_path")
        end
        println("Loading cached D1Q3 norm data from: $h5_path")
        results = load_results_h5(h5_path)
        write_summary(summary_path, results)
        if plot
            plot_results(results, output_pdf)
        end
        println("Summary written to: $summary_path")
        return results
    end

    println("Running D1Q3 full-generator norm sweep for ngrid = $(resolved_ngrid_sweep)")
    println("D1Q3 sweep mode = $sweep_mode")
    println("Output HDF5: $h5_path")
    if plot
        println("Output PDF: $output_pdf")
    end
    println("Summary text: $summary_path")
    println()

    results = NamedTuple[]

    for ngrid_val in resolved_ngrid_sweep
        GC.gc()
        println("="^80)
        @printf("Assembling C_full for ngrid = %d\n", ngrid_val)
        build_elapsed = NaN
        inf_elapsed = NaN
        spectral_elapsed = NaN
        dim = 0
        nnz_val = 0
        spectral_est = NaN
        converged = false
        iterations = 0
        inf_norm_val = NaN
        status = "ok"

        try
            t_build = time()
            C_full, _, core_cfg, case_cfg = build_full_generator(ngrid_val)
            build_elapsed = time() - t_build

            dim = size(C_full, 1)
            nnz_val = nnz(C_full)
            @printf("  lifted dimension = %d\n", dim)
            @printf("  nnz(C_full)      = %d\n", nnz_val)
            @printf("  assembly time    = %.3f s\n", build_elapsed)

            t_norm = time()
            inf_norm_val = opnorm(C_full, Inf)
            inf_elapsed = time() - t_norm
            @printf("  inf norm         = %.8e\n", inf_norm_val)

            if dim <= MAX_SPECTRAL_DIM
                t_spectral = time()
                spectral_est, converged, iterations = power_iteration_spectral_norm(C_full)
                spectral_elapsed = time() - t_spectral
                @printf("  spectral est     = %.8e\n", spectral_est)
                @printf("  spectral time    = %.3f s (converged=%s, iterations=%d)\n", spectral_elapsed, string(converged), iterations)
            else
                status = "spectral estimate skipped (matrix too large for iterative 2-norm estimate in this sweep)"
                println("  spectral est     = skipped (matrix too large for iterative 2-norm estimate in this sweep)")
            end
            @printf("  inf-norm time    = %.3f s\n", inf_elapsed)

            push!(results, (
                ngrid=ngrid_val,
                dim=dim,
                nnz=nnz_val,
                spectral_est=spectral_est,
                inf_norm=inf_norm_val,
                converged=converged,
                iterations=iterations,
                status=status,
                assembly_time_s=build_elapsed,
                inf_norm_time_s=inf_elapsed,
                spectral_time_s=spectral_elapsed,
                Q=core_cfg.Q,
                truncation_order=case_cfg.truncation_order,
            ))
        catch err
            status = sprint(showerror, err)
            println("  ERROR: $status")
            push!(results, (
                ngrid=ngrid_val,
                dim=dim,
                nnz=nnz_val,
                spectral_est=spectral_est,
                inf_norm=inf_norm_val,
                converged=converged,
                iterations=iterations,
                status=status,
                assembly_time_s=build_elapsed,
                inf_norm_time_s=inf_elapsed,
                spectral_time_s=spectral_elapsed,
                Q=3,
                truncation_order=3,
            ))
        end

        write_summary(summary_path, results)
        save_results_h5(h5_path, results)
    end

    println("\n" * "-"^80)
    println(@sprintf("%6s %12s %14s %16s %12s %9s %10s  %s", "ngrid", "dim", "nnz", "spectral_est", "inf_norm", "conv", "iters", "status"))
    for r in results
        println(@sprintf(
            "%6d %12d %14d %16s %12s %9s %10d  %s",
            r.ngrid,
            r.dim,
            r.nnz,
            isfinite(r.spectral_est) ? @sprintf("%.8e", r.spectral_est) : "NaN",
            isfinite(r.inf_norm) ? @sprintf("%.8e", r.inf_norm) : "NaN",
            string(r.converged),
            r.iterations,
            r.status,
        ))
    end

    if plot
        plot_results(results, output_pdf)
    end
    println("Summary written to: $summary_path")
    return results
end

CFSU.run_main_if_script(@__FILE__, main)