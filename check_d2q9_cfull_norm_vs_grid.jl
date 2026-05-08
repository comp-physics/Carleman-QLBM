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

include("src/CLBE/clbe_tg2d_run.jl")
include("src/CLBE/cfull_norm_sweep_utils.jl")

const CFSU = CFullNormSweepUtils

const MAX_SPECTRAL_DIM = 2_000_000
const MAX_ASSEMBLY_DIM = 5_000_000
const DEFAULT_SWEEP_MODE = :arithmetic
const DEFAULT_START_NX = 3
const DEFAULT_STEP_NX = 3
const DEFAULT_END_NX = 6
const DEFAULT_INCREASE_FACTOR = 2
const DEFAULT_NUM_NX = 4

function generate_square_grid_sweep(; start_nx::Int=DEFAULT_START_NX,
    step_nx::Int=DEFAULT_STEP_NX,
    end_nx::Int=DEFAULT_END_NX,
    increase_factor::Int=DEFAULT_INCREASE_FACTOR,
    num_nx::Int=DEFAULT_NUM_NX,
    sweep_mode::Symbol=DEFAULT_SWEEP_MODE)

    if sweep_mode == :arithmetic
        if step_nx <= 0
            error("step_nx must be positive, got $step_nx")
        end
        if end_nx < start_nx
            error("end_nx must be >= start_nx, got start_nx=$start_nx and end_nx=$end_nx")
        end
        return [(nx, nx) for nx in start_nx:step_nx:end_nx]
    elseif sweep_mode == :geometric
        if increase_factor < 2
            error("increase_factor must be >= 2 for geometric sweeps, got $increase_factor")
        end
        if num_nx <= 0
            error("num_nx must be positive, got $num_nx")
        end
        nx_vals = [start_nx * increase_factor^(i - 1) for i in 1:num_nx]
        return [(nx, nx) for nx in nx_vals]
    else
        error("Unsupported sweep_mode=$sweep_mode. Use :arithmetic or :geometric")
    end
end

function tex_output_dir()
    return dirname(CFSU.figure_output_path("d2q9_cfull_norm_vs_grid.pdf"))
end

function summary_output_path()
    return CFSU.data_output_path("d2q9_cfull_norm_vs_grid_summary.txt")
end

function h5_output_path()
    return CFSU.data_output_path("d2q9_cfull_norm_vs_grid.h5")
end

function plot_output_path()
    return CFSU.figure_output_path("d2q9_cfull_norm_vs_grid.pdf")
end

lifted_dimension(Q, truncation_order, nspatial) = sum((Q * nspatial)^k for k in 1:truncation_order)

function power_iteration_spectral_norm(A::SparseMatrixCSC; maxiter::Int=60, tol::Float64=1e-8, seed::Int=1234)
    return CFSU.power_iteration_spectral_norm(A; maxiter=maxiter, tol=tol, seed=seed)
end

function build_full_generator(nx, ny; truncation_order_val=3, coeff_method=:numerical, rho_value=1.0001, boundary_setup=false)
    runtime = prepare_d2q9_carleman_runtime(
        nx=nx,
        ny=ny,
        rho_value=rho_value,
        coeff_method=coeff_method,
        local_truncation_order=truncation_order_val,
        boundary_setup=boundary_setup,
    )

    C_full, bt, _ = build_full_clbe_generator_sparse(
        runtime.setup.symbolic_collision,
        runtime.setup.symbolic_state,
        tau_value,
        Q,
        truncation_order_val,
        poly_order,
        force_factor,
        runtime.w_value,
        runtime.e_value;
        S_lbm=runtime.S_lbm,
        nspatial=runtime.ngrid,
    )

    dropzeros!(C_full)
    return C_full, bt, runtime
end

function write_summary(path, results)
    headers = ["nx", "ny", "dim", "nnz", "spectral_est", "inf_norm", "conv", "iters", "status"]
    rows = [[r.nx, r.ny, r.dim, r.nnz, r.spectral_est, r.inf_norm, string(r.converged), r.iterations, r.status] for r in results]
    CFSU.write_summary_table(path;
        title="D2Q9 full CLBE generator norm sweep",
        intro_lines=[
            "Target matrix: C_full = C - S",
            "Periodic D2Q9 setup, k = 3, numerical coefficients",
        ],
        headers=headers,
        rows=rows,
        left_align=[false, false, false, false, false, false, true, false, true],
    )
end

function save_results_h5(path, results)
    metadata = Dict(
        "target_matrix" => "C_full = C - S",
        "Q" => Int(Q),
        "poly_order" => Int(poly_order),
        "truncation_order" => 3,
        "max_spectral_dim" => Int(MAX_SPECTRAL_DIM),
        "max_assembly_dim" => Int(MAX_ASSEMBLY_DIM),
        "sweep_kind" => "square_grid",
    )
    CFSU.write_results_h5(path, results; metadata=metadata)
end

function load_results_h5(path)
    return CFSU.load_results_h5(path)
end

function plot_results(results, output_pdf)
    CFSU.plot_norm_results(results, output_pdf;
        x_values=[r.nx * r.ny for r in results],
        x_labels=["$(r.nx)×$(r.ny)" for r in results],
        xlabel_text="grid size",
        title_text="D2Q9 full CLBE generator norm vs grid",
        include_inf=true,
        include_spectral=true,
    )
end

function plot_results_from_h5(; h5_path=h5_output_path(), output_pdf=plot_output_path())
    results = load_results_h5(h5_path)
    plot_results(results, output_pdf)
    return results
end

function main(; start_nx::Int=DEFAULT_START_NX,
    step_nx::Int=DEFAULT_STEP_NX,
    end_nx::Int=DEFAULT_END_NX,
    increase_factor::Int=DEFAULT_INCREASE_FACTOR,
    num_nx::Int=DEFAULT_NUM_NX,
    sweep_mode::Symbol=DEFAULT_SWEEP_MODE,
    plot::Bool=(get(ENV, "D2Q9_CFULL_MAKE_PLOT", "1") == "1"),
    load_h5_only::Bool=(get(ENV, "D2Q9_CFULL_LOAD_H5_ONLY", "0") == "1"))

    grid_sweep = generate_square_grid_sweep(
        start_nx=start_nx,
        step_nx=step_nx,
        end_nx=end_nx,
        increase_factor=increase_factor,
        num_nx=num_nx,
        sweep_mode=sweep_mode,
    )
    h5_path = h5_output_path()
    summary_path = summary_output_path()
    output_pdf = plot_output_path()

    if load_h5_only
        if !isfile(h5_path)
            error("Cannot load cached D2Q9 norm data because HDF5 file does not exist: $h5_path")
        end
        println("Loading cached D2Q9 norm data from: $h5_path")
        results = load_results_h5(h5_path)
        write_summary(summary_path, results)
        if plot
            plot_results(results, output_pdf)
        end
        println("Summary written to: $summary_path")
        return results
    end

    println("Running D2Q9 full-generator norm sweep for grids = $(grid_sweep)")
    println("Square-grid sweep mode = $sweep_mode")
    println("Square-grid nx sweep = $([nx for (nx, _) in grid_sweep])")
    println("Output HDF5: $h5_path")
    println("Summary text: $summary_path")
    if plot
        println("Output PDF: $output_pdf")
    end
    println()

    results = NamedTuple[]

    for (nx, ny) in grid_sweep
        GC.gc()
        println("="^80)
        @printf("Assembling D2Q9 C_full for nx = %d, ny = %d\n", nx, ny)
        nspatial = nx * ny
        predicted_dim = lifted_dimension(9, 3, nspatial)
        @printf("  predicted lifted dimension = %d\n", predicted_dim)

        warnings = String[]
        if predicted_dim > MAX_ASSEMBLY_DIM
            warning_msg = @sprintf(
                "WARNING: predicted dim %d exceeds MAX_ASSEMBLY_DIM %d; continuing because this sweep is configured for large-memory / supercomputer runs.",
                predicted_dim,
                MAX_ASSEMBLY_DIM,
            )
            println("  " * warning_msg)
            push!(warnings, warning_msg)
        end

        status = "ok"
        dim = predicted_dim
        nnz_val = 0
        spectral_est = NaN
        converged = false
        iterations = 0
        inf_norm_val = NaN
        build_elapsed = NaN
        inf_elapsed = NaN
        spectral_elapsed = NaN

        try
            t_build = time()
            C_full, _, runtime = build_full_generator(nx, ny)
            build_elapsed = time() - t_build
            dim = size(C_full, 1)
            nnz_val = nnz(C_full)
            @printf("  lifted dimension = %d\n", dim)
            @printf("  nnz(C_full)      = %d\n", nnz_val)
            @printf("  assembly time    = %.3f s\n", build_elapsed)

            t_inf = time()
            inf_norm_val = opnorm(C_full, Inf)
            inf_elapsed = time() - t_inf
            @printf("  inf norm         = %.8e\n", inf_norm_val)
            @printf("  inf-norm time    = %.3f s\n", inf_elapsed)

            if dim > MAX_SPECTRAL_DIM
                warning_msg = @sprintf(
                    "WARNING: lifted dim %d exceeds MAX_SPECTRAL_DIM %d; continuing spectral estimation because this sweep is configured for large-memory / supercomputer runs.",
                    dim,
                    MAX_SPECTRAL_DIM,
                )
                println("  " * warning_msg)
                push!(warnings, warning_msg)
            end

            t_spectral = time()
            spectral_est, converged, iterations = power_iteration_spectral_norm(C_full)
            spectral_elapsed = time() - t_spectral
            @printf("  spectral est     = %.8e\n", spectral_est)
            @printf("  spectral time    = %.3f s (converged=%s, iterations=%d)\n", spectral_elapsed, string(converged), iterations)

            if !isempty(warnings)
                status = join(warnings, " | ")
            end

            results_entry = (
                nx=nx,
                ny=ny,
                predicted_dim=predicted_dim,
                dim=dim,
                nnz=nnz_val,
                spectral_est=spectral_est,
                inf_norm=inf_norm_val,
                converged=converged,
                iterations=iterations,
                status=status,
                coeff_method=runtime.setup.method,
                assembly_time_s=build_elapsed,
                inf_norm_time_s=inf_elapsed,
                spectral_time_s=spectral_elapsed,
            )
            push!(results, results_entry)
        catch err
            status = sprint(showerror, err)
            println("  ERROR: $status")
            push!(results, (
                nx=nx,
                ny=ny,
                predicted_dim=predicted_dim,
                dim=dim,
                nnz=nnz_val,
                spectral_est=spectral_est,
                inf_norm=inf_norm_val,
                converged=converged,
                iterations=iterations,
                status=status,
                coeff_method=:numerical,
                assembly_time_s=build_elapsed,
                inf_norm_time_s=inf_elapsed,
                spectral_time_s=spectral_elapsed,
            ))
        end

        write_summary(summary_path, results)
        save_results_h5(h5_path, results)
    end

    println("\n" * "-"^80)
    for r in results
        println(r)
    end
    if plot
        plot_results(results, output_pdf)
    end
    println("Summary written to: $summary_path")
    return results
end

CFSU.run_main_if_script(@__FILE__, main)