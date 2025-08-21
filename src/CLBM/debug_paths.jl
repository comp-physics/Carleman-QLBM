# Path debugging script for CI
println("🔍 Path Debugging Information")
println("=" ^ 50)

# Environment variables
println("Environment Variables:")
println("  QCFD_HOME = $(get(ENV, "QCFD_HOME", "NOT SET"))")
println("  QCFD_SRC = $(get(ENV, "QCFD_SRC", "NOT SET"))")
println("  PWD = $(pwd())")

# Check if key files exist
println("\nFile Existence Check:")
QCFD_HOME = get(ENV, "QCFD_HOME", "")
QCFD_SRC = get(ENV, "QCFD_SRC", "")

key_files = [
    QCFD_HOME * "/julia_lib/matrix_kit.jl",
    QCFD_HOME * "/visualization/plot_kit.jl", 
    QCFD_SRC * "CLBM/clbm_config.jl",
    QCFD_SRC * "CLBM/timeMarching.jl",
    QCFD_SRC * "LBM/lbm_cons.jl"
]

for file in key_files
    exists = isfile(file)
    status = exists ? "✅" : "❌"
    println("  $status $file")
end

# Directory structure
println("\nDirectory Structure:")
if !isempty(QCFD_HOME) && isdir(QCFD_HOME)
    println("  QCFD_HOME contents:")
    for item in readdir(QCFD_HOME)
        println("    - $item")
    end
end

if !isempty(QCFD_SRC) && isdir(QCFD_SRC)
    println("  QCFD_SRC contents:")
    for item in readdir(QCFD_SRC)
        println("    - $item")
    end
end

println("\n✅ Path debugging complete!")
