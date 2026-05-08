# Validation script to verify CI path fixes
println("🔧 Path Fix Validation")
println("=" ^ 50)

# Test environment variables
println("Environment Variables:")
QCFD_HOME = get(ENV, "QCFD_HOME", "")
QCFD_SRC = get(ENV, "QCFD_SRC", "")
println("  QCFD_HOME = '$QCFD_HOME'")
println("  QCFD_SRC = '$QCFD_SRC'")

# Validate paths are not empty
if isempty(QCFD_HOME)
    error("❌ QCFD_HOME is not set!")
end
if isempty(QCFD_SRC)
    error("❌ QCFD_SRC is not set!")
end

# Test that key files exist
test_files = [
    QCFD_HOME * "/julia_lib/matrix_kit.jl",
    QCFD_HOME * "/visualization/plot_kit.jl",
    QCFD_SRC * "CLBE/clbm_config.jl",
    QCFD_SRC * "CLBE/timeMarching.jl",
    QCFD_SRC * "LBM/lbm_cons.jl"
]

println("\nFile Existence Tests:")
global all_exist = true
for file in test_files
    exists = isfile(file)
    status = exists ? "✅" : "❌"
    println("  $status $file")
    global all_exist = all_exist && exists
end

if all_exist
    println("\n✅ All path validations passed!")
    println("✅ CI environment setup is correct!")
else
    error("❌ Some required files are missing - path setup is incorrect!")
end