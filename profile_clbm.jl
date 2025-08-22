# Profiling script for CLBM code to identify performance bottlenecks

using Profile
using ProfileView

# Load the CLBM code
include("src/CLBM/clbm_config.jl")

using SymPy
using LinearAlgebra
using SparseArrays

include("src/CLBM/coeffs_poly.jl")
include("src/CLBM/collision_sym.jl")
include("src/CLBM/carleman_transferA.jl")
include("src/CLBM/carleman_transferA_ngrid.jl")
include("src/CLBM/LBM_const_subs.jl")
include("src/LBM/lbm_cons.jl")
include("src/LBM/lbm_const_sym.jl")
include("src/LBM/forcing.jl")
include("src/LBM/f_initial.jl")
include("src/CLBM/timeMarching.jl")
include("src/CLBM/CLBM_collision_test.jl")

println("=== CLBM PERFORMANCE PROFILING ===")
println("Configuration: ngrid=$ngrid, use_sparse=$use_sparse")

# Clear any previous profiling data
Profile.clear()

# Run a warmup to compile everything first
println("\n1. Warmup run (compiling functions)...")
@time begin
    w, e, w_val, e_val = lbm_const_sym()
    global w_value = w_val
    global e_value = e_val
    f, omega, u, rho = collision(Q, D, w, e, rho0, lTaylor, lorder2)
    global F1_ngrid, F2_ngrid, F3_ngrid = get_coeff_LBM_Fi_ngrid(poly_order, Q, f, omega, tau_value, ngrid)
    
    if use_sparse
        fT, VT_f, VT = CLBM_collision_test_sparse(Q, omega, f, truncation_order, dt, tau_value, e_value, 5, false)  # Short run
    else
        C, bt, F0 = carleman_C(Q, truncation_order, poly_order, f, omega, tau_value, force_factor, w_value, e_value)
        fT, VT_f, VT = CLBM_collision_test(Q, omega, f, C, truncation_order, dt, tau_value, e_value, 5, false)  # Short run
    end
end

# Now profile the actual run
println("\n2. Profiled run...")

@profile begin
    # Reset globals for clean profiling
    w, e, w_val, e_val = lbm_const_sym()
    global w_value = w_val
    global e_value = e_val
    f, omega, u, rho = collision(Q, D, w, e, rho0, lTaylor, lorder2)
    global F1_ngrid, F2_ngrid, F3_ngrid = get_coeff_LBM_Fi_ngrid(poly_order, Q, f, omega, tau_value, ngrid)
    
    if use_sparse
        println("Profiling SPARSE implementation...")
        fT, VT_f, VT = CLBM_collision_test_sparse(Q, omega, f, truncation_order, dt, tau_value, e_value, n_time, false)
    else
        println("Profiling DENSE implementation...")
        C, bt, F0 = carleman_C(Q, truncation_order, poly_order, f, omega, tau_value, force_factor, w_value, e_value)
        fT, VT_f, VT = CLBM_collision_test(Q, omega, f, C, truncation_order, dt, tau_value, e_value, n_time, false)
    end
end

println("\n=== PROFILING RESULTS ===")

# Print profiling results
Profile.print(format=:flat, sortedby=:count, noisefloor=2.0)

println("\n=== TOP FUNCTIONS BY TIME ===")
Profile.print(format=:flat, sortedby=:count, maxdepth=15)

println("\n=== DETAILED BREAKDOWN ===")

# Get more detailed analysis
data = Profile.fetch()
if !isempty(data)
    println("Total samples: $(length(data))")
    
    # Try to get function-level breakdown
    try
        # Count samples by function name
        func_counts = Dict{String, Int}()
        file_counts = Dict{String, Int}()
        
        for frame in data
            if frame != 0
                # Get the stack frame info
                frame_info = Base.StackTraces.lookup(frame)
                if !isempty(frame_info)
                    sf = frame_info[1]
                    func_name = string(sf.func)
                    file_name = string(sf.file)
                    
                    # Count by function
                    if haskey(func_counts, func_name)
                        func_counts[func_name] += 1
                    else
                        func_counts[func_name] = 1
                    end
                    
                    # Count by file
                    if haskey(file_counts, file_name)
                        file_counts[file_name] += 1
                    else
                        file_counts[file_name] = 1
                    end
                end
            end
        end
        
        # Sort and display top functions
        sorted_funcs = sort(collect(func_counts), by=x->x[2], rev=true)
        println("\nTop functions by sample count:")
        for (i, (func, count)) in enumerate(sorted_funcs[1:min(20, end)])
            pct = round(count / length(data) * 100, digits=1)
            println("  $i. $func: $count samples ($pct%)")
        end
        
        # Sort and display top files
        sorted_files = sort(collect(file_counts), by=x->x[2], rev=true)
        println("\nTop files by sample count:")
        for (i, (file, count)) in enumerate(sorted_files[1:min(15, end)])
            pct = round(count / length(data) * 100, digits=1)
            # Only show relevant files (skip Julia internals)
            if occursin("CLBM", file) || occursin("LBM", file) || occursin("carleman", file) || occursin("timeMarching", file)
                println("  $i. $(basename(file)): $count samples ($pct%)")
            end
        end
        
    catch e
        println("Error in detailed analysis: $e")
    end
else
    println("No profiling data collected. Try running with a longer simulation.")
end

println("\n=== PROFILING COMPLETE ===")
println("To view the interactive profiling results, you can run:")
println("using ProfileView; ProfileView.view()")
println("Or save a flame graph with:")
println("using FlameGraphs, ProfileCanvas; ProfileCanvas.html_file(\"profile.html\", flamegraph())")
