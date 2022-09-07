abstract type ContinuousStateModel <: BranchModel end #Unclear that we need this mid layer.

include("gaussian_partition.jl")
include("brownian_motion.jl")
include("ornstein_uhlenbeck.jl")
