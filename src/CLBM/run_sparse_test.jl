# Simple test runner for sparse vs dense comparison
# Run this with: julia run_sparse_test.jl

println("=" ^ 60)
println("SPARSE VS DENSE CARLEMAN MATRIX TEST")
println("=" ^ 60)

try
    include("test_sparse_vs_dense.jl")
    println("\n" * "=" ^ 60)
    println("TEST COMPLETED SUCCESSFULLY!")
    println("=" ^ 60)
catch e
    println("\n" * "=" ^ 60)
    println("TEST FAILED!")
    println("Error: ", e)
    println("=" ^ 60)
    rethrow(e)
end
