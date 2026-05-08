using Printf
using LinearAlgebra
using SparseArrays
using Dates

if !haskey(ENV, "QCFD_HOME")
    ENV["QCFD_HOME"] = pwd()
end
if !haskey(ENV, "QCFD_SRC")
    ENV["QCFD_SRC"] = joinpath(pwd(), "src") * "/"
end

include("src/CLBE/clbe_multigrid_run.jl")

const VERIFY_NGRID = [1, 3, 6, 12, 24, 48]

function summary_output_path()
    outdir = joinpath(pwd(), "data")
    mkpath(outdir)
    return joinpath(outdir, "d1q3_generator_path_equality_summary.txt")
end

function verify_one_case(ngrid_val; truncation_order_val=3, coeff_method=:numerical, dt_val=1.0)
    core_cfg, case_cfg = build_d1q3_multigrid_configs(
        comparison_ngrid=ngrid_val,
        local_use_sparse=true,
        local_n_time=2,
        local_truncation_order=truncation_order_val,
        coeff_method=coeff_method,
        initial_condition=:legacy,
        dt_override=dt_val,
    )

    runtime = prepare_d1q3_runtime(core_cfg, case_cfg)
    phi_ini = multigrid_initial_condition(ngrid_val; initial_condition=case_cfg.initial_condition, u_ini=case_cfg.u_ini)
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

    V0 = Float64.(carleman_V(phi_ini, case_cfg.truncation_order))
    VT_expected = V0 .+ dt_val .* (C_full * V0 .+ bt)

    phiT, VT = timeMarching_state_CLBM_sparse(
        runtime.omega,
        runtime.f,
        core_cfg.tau_value,
        core_cfg.Q,
        case_cfg.truncation_order,
        dt_val,
        phi_ini,
        2;
        S_lbm=S_lbm,
        nspatial=ngrid_val,
        integrator=:euler,
    )

    vt_diff = VT[:, 2] .- VT_expected
    phi_diff = phiT[:, 2] .- VT_expected[1:length(phi_ini)]

    return (
        ngrid=ngrid_val,
        dim=size(C_full, 1),
        nnz=nnz(C_full),
        vt_max_abs=maximum(abs.(vt_diff)),
        vt_l2=norm(vt_diff),
        phi_max_abs=maximum(abs.(phi_diff)),
        phi_l2=norm(phi_diff),
        F1_sparse=issparse(F1_ngrid),
        F2_sparse=(F2_ngrid === nothing ? true : issparse(F2_ngrid)),
        F3_sparse=(F3_ngrid === nothing ? true : issparse(F3_ngrid)),
    )
end

function write_summary(path, results)
    open(path, "w") do io
        println(io, "D1Q3 generator path equality verification")
        println(io, "Generated: ", Dates.format(now(), dateformat"yyyy-mm-dd HH:MM:SS"))
        println(io, "Checks that the production time-marching path and norm-sweep path use the same shared sparse numerical coefficient-lifting and full-generator assembly.")
        println(io)
        println(io, @sprintf("%6s %12s %12s %14s %14s %14s %14s %8s %8s %8s", "ngrid", "dim", "nnz", "vt_max_abs", "vt_l2", "phi_max_abs", "phi_l2", "F1sp", "F2sp", "F3sp"))
        for r in results
            println(io, @sprintf(
                "%6d %12d %12d %14.6e %14.6e %14.6e %14.6e %8s %8s %8s",
                r.ngrid,
                r.dim,
                r.nnz,
                r.vt_max_abs,
                r.vt_l2,
                r.phi_max_abs,
                r.phi_l2,
                string(r.F1_sparse),
                string(r.F2_sparse),
                string(r.F3_sparse),
            ))
        end
    end
end

function main()
    summary_path = summary_output_path()
    results = NamedTuple[]

    println("Verifying D1Q3 production-path equality for ngrid = $(VERIFY_NGRID)")
    for ngrid_val in VERIFY_NGRID
        println("="^80)
        @printf("Checking ngrid = %d\n", ngrid_val)
        result = verify_one_case(ngrid_val)
        push!(results, result)
        @printf("  dim            = %d\n", result.dim)
        @printf("  nnz            = %d\n", result.nnz)
        @printf("  vt_max_abs     = %.6e\n", result.vt_max_abs)
        @printf("  vt_l2          = %.6e\n", result.vt_l2)
        @printf("  phi_max_abs    = %.6e\n", result.phi_max_abs)
        @printf("  phi_l2         = %.6e\n", result.phi_l2)
        @printf("  sparse F mats  = (%s, %s, %s)\n", string(result.F1_sparse), string(result.F2_sparse), string(result.F3_sparse))
        write_summary(summary_path, results)
    end

    println("\nSummary written to: $summary_path")
end

main()