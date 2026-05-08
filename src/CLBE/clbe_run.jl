println("ℹ️  `src/CLBE/clbe_run.jl` is kept as a compatibility wrapper.")
println("   The primary operational driver is `src/CLBE/clbe_multigrid_run.jl`.")

include("clbe_multigrid_run.jl")
