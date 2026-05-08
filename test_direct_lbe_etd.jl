using LinearAlgebra
using SparseArrays
using Printf

const REPO_ROOT = @__DIR__
const SRC_ROOT = joinpath(REPO_ROOT, "src")

ENV["QCFD_HOME"] = get(ENV, "QCFD_HOME", REPO_ROOT)
ENV["QCFD_SRC"] = get(ENV, "QCFD_SRC", SRC_ROOT * "/")

module DirectLBETestModule
include(joinpath(Main.SRC_ROOT, "LBE", "direct_LBE.jl"))
end

module D1Q3ETDTestModule
include(joinpath(Main.SRC_ROOT, "CLBE", "clbe_multigrid_run.jl"))
end

module D2Q9ETDTestModule
include(joinpath(Main.SRC_ROOT, "CLBE", "clbe_tg2d_run.jl"))
end

function passfail(ok::Bool)
    return ok ? "PASS" : "FAIL"
end

function print_banner(title)
    println()
    println("="^72)
    println(title)
    println("="^72)
end

function run_linear_etd_accuracy_test()
    print_banner("Linear ETD accuracy test")

    phi0 = [1.0, -2.0]
    S_lbm = spzeros(2, 2)
    F1_ngrid = sparse(Diagonal([-1.0, 0.5]))
    F2_ngrid = nothing
    F3_ngrid = nothing
    dt = 0.3

    history = DirectLBETestModule.timeMarching_direct_LBE_ngrid(
        phi0,
        dt,
        2,
        F1_ngrid,
        F2_ngrid,
        F3_ngrid;
        S_lbm=S_lbm,
        integrator=:exponential_euler,
    )
    exact = exp(dt * Matrix(F1_ngrid)) * phi0
    step_err = norm(history[:, 2] - exact)
    ok = isfinite(step_err) && step_err <= 1e-12

    println("ETD step error = ", step_err)
    println("$(passfail(ok)) linear_etd_accuracy")
    return ok
end

function run_d1q3_etd_smoke_test()
    print_banner("D1Q3 ETD smoke test")

    phiT_lbe, phiT_clbm, VT, avg_abs_err, avg_rel_err = D1Q3ETDTestModule.main(
        comparison_ngrid=6,
        local_n_time=40,
        local_truncation_order=3,
        l_plot=false,
        coeff_method=:numerical,
        initial_condition=:legacy,
        dt_override=1.0,
        integrator=:matrix_exponential,
        direct_lbe_integrator=:exponential_euler,
    )

    lbe_finite = all(isfinite, phiT_lbe)
    clbe_finite = all(isfinite, phiT_clbm)
    max_abs_err = maximum(avg_abs_err)
    ok = lbe_finite && clbe_finite && isfinite(max_abs_err)

    println("direct LBE finite = ", lbe_finite)
    println("CLBE finite = ", clbe_finite)
    println("max domain-averaged absolute error = ", max_abs_err)
    println("$(passfail(ok)) d1q3_etd_smoke")
    return ok
end

function run_d2q9_etd_smoke_test()
    print_banner("D2Q9 TG ETD smoke test")

    phiT_ref, phiT_clbm, VT, dist_abs_err, density_error_norm, vel_abs_err = D2Q9ETDTestModule.main(
        nx=3,
        ny=3,
        amplitude=0.02,
        local_n_time=10,
        l_plot=false,
        boundary_setup=false,
        coeff_method=:numerical,
        local_truncation_order=3,
        reference_model=:direct_lbe,
        integrator=:matrix_exponential,
        direct_lbe_integrator=:exponential_euler,
    )

    ref_finite = all(isfinite, phiT_ref)
    clbe_finite = all(isfinite, phiT_clbm)
    max_dist_err = maximum(dist_abs_err)
    max_density_err = maximum(density_error_norm)
    max_velocity_err = maximum(vel_abs_err)
    ok = ref_finite && clbe_finite && isfinite(max_dist_err) && isfinite(max_density_err) && isfinite(max_velocity_err)

    println("direct LBE finite = ", ref_finite)
    println("CLBE finite = ", clbe_finite)
    println("max distribution absolute error = ", max_dist_err)
    println("max density error norm = ", max_density_err)
    println("max velocity error norm = ", max_velocity_err)
    println("$(passfail(ok)) d2q9_etd_smoke")
    return ok
end

function main()
    results = Dict(
        :linear_etd_accuracy => run_linear_etd_accuracy_test(),
        :d1q3_etd_smoke => run_d1q3_etd_smoke_test(),
        :d2q9_etd_smoke => run_d2q9_etd_smoke_test(),
    )

    print_banner("ETD test summary")
    all_ok = true
    for key in (:linear_etd_accuracy, :d1q3_etd_smoke, :d2q9_etd_smoke)
        ok = results[key]
        all_ok &= ok
        println(rpad(String(key), 28), " : ", passfail(ok))
    end

    if !all_ok
        error("At least one ETD test failed.")
    end

    println("All ETD tests passed.")
end

main()