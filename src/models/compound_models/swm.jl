#BM: Most importantly: This needs to be tested on simple numerical examples.
#Test 1: set up a tree with a single partition, run felsenstein for N different models, and manually compute what the mixture likelihood should be. Compare to SWMModel LL.
#Test 2: Check that the LLs, weighting, etc, at all nodes are equal after the down pass (eg. when computing node marginals).

#Slight annoyance with the current setup: The Partition needs to know the weights, because site_LLs() needs to know its weights, and that doesn't see the model.
#But if we change the weights, we kinda want them changed everywhere, which we'd need to do manually.
#Solution for now: have the model and the partition both store the weights, and backward! and forward! will update the Partition from the model.
#This doesn't work, because the total LL is called on the root partition, and forward! and backward! don't touch this.
#Should they change source.weights as well? Ok, lets try that. This requires forcing a prop for the root partition, which is a bit annoying.

#Another way of thinking about this is along the same lines as the root freqs. Leave them out of the model, and force you to set them manually at the root partitions.
#Need to think about how this will interact with the blopt, NNI, ancestors, etc. Those nodes don't know their weights, but do they get that info from the down pass?
#Maybe the weights from above just tweak the site scalings, and then these get used in all site_LL calcs to decide the weight of contribution from each partition
#for each site?

export SWMModel
mutable struct SWMModel <: BranchModel
    models::Vector{<:BranchModel}
    weights::Vector{Float64}
    function SWMModel(models::Vector{<:BranchModel})
        new(models, ones(length(models)) ./ length(models))
    end
    #Convenience constructor: takes a single model and an array of rates
    function SWMModel(model::M, rs::Vector{Float64}) where {M <: BranchModel}
        models = (vcat(model, [deepcopy(model) for x in rs[2:end]]))
        for (i,m) in enumerate(models)
            m.r = rs[i]
        end
        new(models, ones(length(models)) ./ length(models))
    end
end


export SWMPartition
#Parts is a vector containing the various components
#Weights is a num_components by num_sites matrix
mutable struct SWMPartition{PType} <: Partition where {PType <: MultiSitePartition}
    parts::Vector{PType}
    weights::Vector{Float64}
    sites::Int #Assumption: sites is the length of the vector returned by site_LLs(part). This is where the mixing happens. ie. sites are the things that are independent under the SWM
    states::Int
    models::Int
    #constructor where sizes are sufficient to initialize parts
    function SWMPartition{PType}(parts::Vector{PType}) where {PType <: MultiSitePartition}
        l = length(parts)
        new{PType}(parts, ones(l)./l, parts[1].sites, parts[1].states, l)
    end
    function SWMPartition{PType}(part::PType, n_parts::Int) where {PType <: MultiSitePartition}
        SWMPartition{PType}([copy_partition(part) for i in 1:n_parts])
    end
    function SWMPartition{PType}(parts::Vector{PType}, weights::Vector{Float64}, sites::Int, states::Int, models::Int) where {PType <: MultiSitePartition}
        @assert length(parts) == length(weights)
        new{PType}(parts, weights, sites, states, models)
    end
end

#Overloading the copy_partition to avoid deepcopy.
function copy_partition(src::SWMPartition{PType}) where {PType <: MultiSitePartition}
    return SWMPartition{PType}(copy_partition.(src.parts), copy(src.weights), src.sites, src.states, src.models)
end

function combine!(dest::SWMPartition{PType},src::SWMPartition{PType}) where {PType<:MultiSitePartition}
    for i in 1:length(dest.parts)
        combine!(dest.parts[i], src.parts[i])
    end
end

function identity!(dest::SWMPartition)
    for part in dest.parts
        identity!(part)
    end
end

#"backward" refers to propagating up the tree: ie time running backwards.
function backward!(dest::SWMPartition{PType},
        source::SWMPartition{PType},
        model::SWMModel,
        node::FelNode) where {PType<:MultiSitePartition}
    dest.weights .= model.weights
    source.weights .= model.weights
    for i in 1:length(dest.parts)
        backward!(dest.parts[i], source.parts[i], model.models[i], node)
    end
end
#"forward" refers to propagating up the tree: ie time running forwards.
function forward!(dest::SWMPartition{PType},
        source::SWMPartition{PType},
        model::SWMModel,
        node::FelNode) where {PType<:MultiSitePartition}
    source.weights .= model.weights
    dest.weights .= model.weights
    for i in 1:length(dest.parts)
        forward!(dest.parts[i], source.parts[i], model.models[i], node)
    end
end


function eq_freq(model::SWMPartition)
    error("eq_freq() not yet implemented for $(typeof(model)). Required for automatically setting root parent messages.")
    #return Vector{T}
end

#-----BM: Everything below here needs checking and testing-----
#Note: much of the code for sampling is *not* going to be very performant, because it propogates over
#partitions that we know, at the root, are already irrelevant. Most applications don't need this to be very fast,
#and it might be easier to make a different `SWMSimPartition` for efficient sampling, but not inference?

function weighted_prob_grid(part::SWMPartition{PType}) where {PType <: MultiSitePartition}
    LL_grid = zeros(part.models, part.sites)
    for i in 1:part.models
        LL_grid[i,:] .= site_LLs(part.parts[i]) .+ log(part.weights[i])
    end
    #Maybe this should be sums...
    maxis = maximum(LL_grid, dims = 1)[:]
    LL_grid .-= maxis'
    LL_grid .= exp.(LL_grid)
    return LL_grid, maxis
end

export SWM_prob_grid
#Very duplicatey. Could replace the above (which is used in a few core calculations) with this, and then taking on the weights after.
"""
    SWM_prob_grid(part::SWMPartition{PType}) where {PType <: MultiSitePartition}

Returns a matrix of probabilities for each site, for each model (in the probability domain - not logged!) as well as the log probability offsets
"""
function SWM_prob_grid(part::SWMPartition{PType}) where {PType <: MultiSitePartition}
    LL_grid = zeros(part.models, part.sites)
    for i in 1:part.models
        LL_grid[i,:] .= site_LLs(part.parts[i])
    end
    maxis = maximum(LL_grid, dims = 1)[:]
    LL_grid .-= maxis'
    LL_grid .= exp.(LL_grid)
    return LL_grid, maxis
end

#This calculates the weights for each component, given the log likelihoods of each component at each site.
#And then samples which partition to use from these.
#And sets the scaling consts for each partition to 0 if the site was picked, and -Inf if not picked.
#These will get propogated by forward!
#Should give the right behavior for sample_down! if the root partition is set up correctly,
#But also for endpoint_conditioned stuff. Need to check those dictionaries...
function sample_partition!(partition::SWMPartition{PType}) where {PType <: MultiSitePartition}
    LL_grid,_ = weighted_prob_grid(partition)
    for i = 1:partition.sites
        LL_grid[:, i] .= one_hot_sample(LL_grid[:, i])
    end
    for i in 1:partition.models
        partition.parts[i].scaling .= log.(LL_grid[i,:])
    end
    for p in partition.parts
        sample_partition!(p)
    end
end

#We need to formalize a weighted mixture to get anything out of these.
#Will be useful for marginal dicts, as well partition2obs, and possibly sampling?
#Clear for DiscretePartitions. Messy but defined for eg. Gaussians (we could make it return a mixture from Distributions)
"""
    mix(swm_part::SWMPartition{PType} ) where {PType <: MultiSitePartition}

`mix` collapses a Site-Wise Mixture partition to a single component partition, weighted by the site-wise likelihoods for each component, and the init weights.
Specifically, it takes a `SWMPartition{Ptype}` and returns a `PType`.
You'll need to have this implemented for certain helper functionality if you're playing with new kinds of SWMPartitions that aren't mixtures of `DiscretePartitions`.
"""
function mix(swm_part::SWMPartition{PType} ) where {PType <: DiscretePartition}
    prob_grid,_ = weighted_prob_grid(swm_part)
    out = copy_partition(swm_part.parts[1])
    out.scaling .= 0.0
    out.state .= 0.0
    for i in 1:swm_part.models
       for j in 1:swm_part.sites
            v = swm_part.parts[i].state[:,j]
            v ./= sum(v)
           out.state[:,j] += prob_grid[i,j] .* v
       end
    end
    return out
end

function mix(swm_part::SWMPartition{PType} ) where {PType <: MultiSitePartition}
    @error "Implement mix() for $(typeof(swm_part)) to be able to do this."
end

function partition2obs(part::SWMPartition{PType}) where {PType <: MultiSitePartition}
    partition2obs(mix(part))
end

function obs2partition!(dest::SWMPartition{PType},seq::String) where {PType <: MultiSitePartition}
    for x in dest.parts
        obs2partition!(x, seq)
    end
end

function site_LLs(part::SWMPartition)
    LL_grid,maxis = weighted_prob_grid(part)
    return log.(sum(LL_grid, dims = 1)[:]) .+ maxis
end

function eq_freq_from_template(model::SWMModel,partition_template::SWMPartition{PType}) where {PType<:MultiSitePartition}
    out_partition = copy_partition(partition_template)
    for i in 1:length(out_partition.parts)
        out_partition.parts[i] = eq_freq_from_template(model.models[i], out_partition.parts[i])
    end
    return out_partition
end

#=
#Need to implement max_partition! for SWMPartitions.
#We'll first take the max over which category is used, and then, conditioned on that, we take the max over the partition.
=#
function max_partition!(partition::SWMPartition{PType}) where {PType <: MultiSitePartition}
    LL_grid,_ = weighted_prob_grid(partition)
    for i = 1:partition.sites
        m = argmax(LL_grid[:, i])
        LL_grid[:, i] .= 0.0
        LL_grid[m, i] = 1.0
    end
    for i in 1:partition.models
        partition.parts[i].scaling .= log.(LL_grid[i,:])
    end
    for p in partition.parts
        max_partition!(p)
    end
end

#Reflection: It kinda seems like we need a "reweight" function that sets SWM-wide site scaling constants to the max of the scalings for each partition, per-site,
#sets all the individual partition scalings to zero, and reweights the states for each site accordingly. Then state/category distributions can just be read off the
#current state values, and the site_LLs function becomes trivial.
#We'll need to handle the SWM-wide scalings during fwd,bck,combine etc, but that should be easy?
#Unclear is this would give speedups, but I think it would reduce bugs. Might also reduce allocs?
#We could also just set the sub-partitions to be equal to the max over all of them, and not need to store SWM-wide scalings?