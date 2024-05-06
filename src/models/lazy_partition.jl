#If you define a new partition type, that isn't a kind of DiscretePartition, then you need to define:
#combine!(dest::NewPartitionType,src::NewPartitionType)
#Combines src and dest and sets into dest
#identity!(dest::NewPartitionType)
#Sets dest to the identity value such that combine!(dest, other) == other
#site_LLs(dest::MyNewPartition)
#Gets a site-wise log-likelihood score

export LazyPartition
mutable struct LazyPartition{PType} <: Partition where {PType <: Partition}
    partition::Union{PType, Nothing}
    memoryblocks::Vector{PType} #A stack that is meant for LazyPartitions in the same tree to share
    static::Bool #Does the LazyPartition keep its partition during a message passing wave? Defaults to false
    obs #Leaf nodes can store their observations here in a preferably compact data structure

    function LazyPartition{PType}(partition) where {PType <: Partition}
        new(partition, Vector{PType}(), false)
    end

    function LazyPartition{PType}(partition, memoryblocks) where {PType <: Partition}
        new(partition, memoryblocks, false)
    end
end

function copy_partition(src::LazyPartition{PType}) where {PType <: Partition}
    if isnothing(src.partition)
        partition = nothing
    elseif isempty(src.memoryblocks)
        partition = copy_partition(src.partition)
    else #Avoid allocating new memory
        partition = pop!(src.memoryblocks)
        copy_partition_to!(partition, src.partition)
    end
    return LazyPartition{PType}(partition, src.memoryblocks) #We want to have a reference to the same stack
end

function combine!(dest::LazyPartition{PType}, src::LazyPartition{PType}) where {PType <: Partition}
    (isnothing(src.partition) || isnothing(dest.partition)) && throw(ArgumentError("The partition field in both the source and destination LazyPartitions must be defined to combine them."))
    combine!(dest.partition, src.partition)
    safe_release_partition!(src)
end


"""
The idea with backward!{LazyPartition} is to avoid having to call the GC too often,
therefore, we need access to shared memoryblocks during the message passing wave. This is only really necessary if PType <: MultisitePartition,
but for consistency, we might as well do it for all types?
"""

"""
I store:
    - References to partitions in a global stack
    - Observations in the obs field of the LazyPartition

To achieve a maximum amount of memoryblocks available during felsenstein!, the workflow of backward! and combine! is this:
    - We pop a memoryblock of the stack and put it in dest.partition prior to a backward!.
    Now, we perform the backward!, then push our source.partition to the memoryblock stack
    - After a combine!, we push src.partition to the memoryblock stack
"""

function backward!(
    dest::LazyPartition{PType},
    source::LazyPartition{PType},
    model::BranchModel,
    node::FelNode,
) where {PType <: Partition}
    if isleafnode(node) && isdefined(source, :obs)
        source.partition = pop!(source.memoryblocks)
        #Transform source.obs to the appropriate format
        lazy_obs2partition!(source.partition, source.obs)
        # In the case of CodonPartition, we need to enforce scaling being zeros (which we do with lazy_obs2partition)
    end
    dest.partition = pop!(source.memoryblocks)
    (isnothing(source.partition) || isnothing(dest.partition)) && throw(ArgumentError("The partition field in the source and dest LazyPartition must be something in order to propagate it backwards."))
    backward!(dest.partition, source.partition, model, node)
    safe_release_partition!(source)
end

function forward!(
    dest::LazyPartition{PType},
    source::LazyPartition{PType},
    model::BranchModel,
    node::FelNode,
) where {PType <: Partition}
    forward!(dest.partition, source.partition, model, node)
    safe_release_partition!(source)
end

function site_LLs(dest::LazyPartition)
    result = site_LLs(dest.partition)
    safe_release_partition!(dest) #All necessary information from dest is extracted into result
    return result
end

function copy_partition_to!(dest::LazyPartition{PType}, src::LazyPartition{PType}) where {PType <: Partition}
    dest.partition = src.partition
    src.partition = nothing #Hmm, we might not always want to wipe the src?
end

"""
Should be run on a tree containing LazyPartitions before running `felsenstein!`. Sorts for a minimal count of active partitions during a felsenstein!
Note: since felsenstein! uses a stack, we want to avoid having long node.children[1].children[1]... chains
"""
function lazysort!(node)
    if isleafnode(node)
        return 1
    end
    maximum_active_partitions = [lazysort!(child) for child in node.children]
    perm = sortperm(maximum_active_partitions)
    node.children = node.children[perm]
    for (i, idx) in enumerate(perm)
        maximum_active_partitions[idx] += length(node.children) - i
    end
    return maximum(maximum_active_partitions)
end

function safe_release_partition!(src::LazyPartition)
    if !src.static
        push!(src.memoryblocks, src.partition)
        src.partition = nothing
    end
end

#Store the obs in the dest LazyPartition
function obs2partition!(dest::LazyPartition, obs)
    dest.obs = obs
end

export lazyprep!
"""
    lazyprep!(tree::FelNode, eq_message::Vector{<:Partition}; partition_list = 1:length(tree.message))

1. Perform a lazysort! on tree to obtain the optimal tree for a lazy felsenstein! prop
2. Fix tree.parent_message to the equilibrium message
3. Add sufficiently many partitions to memoryblocks needed for a felsenstein! prop
"""
function lazyprep!(tree::FelNode, eq_message::Vector{<:Partition}; partition_list = 1:length(tree.message))
    @assert length(eq_message) == length(partition_list)
    maximum_active_partitions = lazysort!(tree)
    for (p, part) in enumerate(partition_list)
        parent_message = tree.parent_message[part]
        parent_message.partition = eq_message[p]
        parent_message.static = true
        for _ = 1:maximum_active_partitions+1
            push!(parent_message.memoryblocks, partition_from_template(eq_message[p]))
        end
    end
    return maximum_active_partitions
end

#Return a new instance of T with the same shape as partition_template. Here, we don't care about the values stored within the partition.
#Fallback method. This should be overloaded with usage of undef to increase performance
function partition_from_template(partition_template::T) where {T <: Partition}
    return copy_partition(partition_template)
end

function partition_from_template(partition_template::T) where {T <: DiscretePartition}
    states, sites = partition_template.states, partition_template.sites
    return T(Array{Float64, 2}(undef, states, sites), states, sites, Array{Float64, 1}(undef, sites))
end

function partition_from_template(partition_template::SWMPartition{PType}) where {PType <: MultiSitePartition}
    return SWMPartition{PType}(partition_from_template.(partition_template.parts), 
        copy(partition_template.weights),
        partition_template.sites,
        partition_template.states,
        partition_template.models)
end