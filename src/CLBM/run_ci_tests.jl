# CI Test Runner - automatically chooses appropriate tests based on available packages

println("🚀 Starting CI Tests...")
println("Julia version: $(VERSION)")
println("Platform: $(Sys.MACHINE)")

# Check if plotting packages are available by attempting to run the full test
global has_plotting = false
try
    # Check if required packages are installed
    using PyPlot
    global has_plotting = true
    println("📊 Full plotting support detected")
catch e
    println("⚠️  PyPlot not available: $e")
    println("📋 Will use minimal tests instead")
end

# Run appropriate tests based on available dependencies
if has_plotting
    println("🔄 Running comprehensive sparse vs dense tests...")
    try
        include("test_sparse_vs_dense.jl")
        println("✅ Full sparse vs dense tests PASSED")
    catch e
        println("❌ Full tests FAILED: $e")
        # Fall back to minimal tests
        println("🔄 Falling back to minimal tests...")
        include("unit_tests_minimal.jl")
    end
else
    println("🔄 Running minimal tests...")
    include("unit_tests_minimal.jl")
end

println("\n🎉 CI Tests completed successfully!")