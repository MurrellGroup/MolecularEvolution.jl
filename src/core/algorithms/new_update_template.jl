#=
A guide for the AbstractUpdate.jl interface.
- See ?MolecularEvolution.AbstractUpdate
- An AbstractUpdate can be really anything, therefore
- We show how you can make flexible use of StandardUpdate
- See ?StandardUpdate
=#

#--- nni_selection::Function ---
#=
The two typical options are
1. Optimization: nni_selection = argmax
2. Metropolis sampling: nni_selection = softmax_sampler
But you can also implement it as
3. Any function that takes a vector of log-likelihoods, and returns an index.
=#

#--- branchlength_modifier::UnivariateModifier ---
#=
The two typical options are
1. Optimization: branchlength_sampler = GoldenSectionOpt() or BrentsMethodOpt()
2. Metropolis sampling: branchlength_sampler = BranchlengthSampler(<symmetric proposal>, <prior>)
3. Making a new UnivariateModifier
=#
struct MyBranchlengthModifier <: UnivariateModifier end

function univariate_modifier(f, modifier::MyBranchlengthModifier; kwargs...)
    error("univariate_modifier() not yet implemented for $(typeof(modifier)). Required for branchlength_update!.")
end

#--- root_update::RootUpdate ---
#=
The two typical options are
1. Optimization: root_update::RootOpt
2. Metropolis sampling: root_update::RootSample
3. Making a new RootUpdate
=#
#--- 1 ---
struct MyRootOpt <: RootOpt end

#either
function Base.length(root_opt::MyRootOpt)
    error("length() not yet implemented for $(typeof(root_opt)). Required for root_opt!.")
end
function (root_opt::MyRootOpt)(message::Vector{<:Partition})
    error()
end

#or
function (root_opt::MyRootOpt)(tree::FelNode, models, partition_list, node_message::Vector{<:Partition}, temp_message::Vector{<:Partition})
    error("")
end

#--- 2 ---
struct MyRootSample <: RootSample end
#either implement the metropolis_step interface and...
Base.length(root_sample::MyRootSample) = error() # the number of consecutive samples of root state and position for a single update call

#or
function (root_sample::MyRootSample)(tree::FelNode, models, partition_list, node_message::Vector{<:Partition}, temp_message::Vector{<:Partition})
    error("")
end

#--- 3 ---
struct MyRootUpdate <: RootUpdate end

function (my_root_update::MyRootUpdate)(tree::FelNode, models, partition_list, node_message::Vector{<:Partition}, temp_message::Vector{<:Partition})
    error("root_update!() not yet implemented for $(typeof(my_root_update)). Required for root_update!.")
end


#--- models_update::ModelsUpdate ---
# See example in docs "Updating a phylogenetic tree"
#TODO: write this up properly

#I also have some convenience structs for the most common cases