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
    temp::Union{PType, Nothing} #This is a temporary memory location for the message passing wave.
    memoryblocks::Vector{PType} #A stack that is meant for LazyPartitions to share

    function LazyPartition{PType}(partition) where {PType <: Partition}
        new(partition, deepcopy(partition)) #Maybe there is a better way to get a new instance of the same type?
    end
end

function combine!(dest::LazyPartition{PType}, src::LazyPartition{PType}) where {PType}
    (isnothing(src.partition) || isnothing(dest.partition)) && throw(ArgumentError("The partition field in both the source and destination LazyPartitions must be defined to combine them."))
    combine!(dest.partition, src.partition)
    push!(src.memoryblocks, src.partition)
    src.partition, src.temp = nothing, nothing
end


"""
The idea with backward!{LazyPartition} is to avoid having to call the GC too often,
therefore, we always need a temp memory location in the message passing wave. This is only really necessary if PType <: MultisitePartition,
but for consistency, we might as well do it for all types?
"""

"""
I store:
    - References to partitions in a global stack
    - Sequences at leafnodes

To achieve a maximum amount of memoryblocks available during felsenstein!, the workflow of backward! and combine! is this:
    - We pop a memoryblock of the stack and put it in source.temp prior to a backward!.
    Now, we perform the backward!, then push our source.partition to the memoryblock stack
    - We let our temps be nothing during the combine! ops as they are irrelevant here.
"""

function backward!(
    dest::LazyPartition{PType},
    source::LazyPartition{PType},
    model::BranchModel,
    node::FelNode,
) where {PType}
    if isleafnode(node)
        source.partition = pop!(source.memoryblocks)
        lazy_obs2partition!(source.partition, node.node_data["obs"])
        # In the case of CodonPartition, we need to enforce scaling being zeros (which we do with lazy_obs2partition)
    end
    source.temp = pop!(source.memoryblocks)
    (isnothing(source.partition) || isnothing(source.temp)) && throw(ArgumentError("The partition and temp fields in the source LazyPartition must be defined to propagate it backwards."))
    backward!(source.temp, source.partition, model, node)
    push!(source.memoryblocks, source.partition)
    dest.partition = source.temp
    source.partition, source.temp = nothing, nothing
end

function site_LLs(dest::LazyPartition)
    return site_LLs(dest.partition)
end

function copy_partition_to!(dest::LazyPartition{PType}, src::LazyPartition{PType}) where {PType}
    dest.partition, dest.temp = src.partition, src.temp
    src.partition, src.temp = nothing, nothing #Hmm, we might not always want to wipe the src?
end

"""
Sorts for a minimal count of active partitions during a felsenstein!
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

#Might not keep this, just for testing the principal
function populate_tree!(
    tree::FelNode,
    starting_message::LazyPartition,
    names,
    data;
    init_all_messages = true,
    tolerate_missing = 1, #0 = error if missing; 1 = warn and set to missing data; 2 = set to missing data
)
    if init_all_messages
        internal_message_init!(tree, starting_message)
    else
        tree.parent_message = deepcopy(starting_message)
    end
    name_dic = Dict(zip(names, 1:length(names)))
    for n in getleaflist(tree)
        if haskey(name_dic, n.name)
            n.node_data = Dict("obs" => data[name_dic[n.name]])
        else
            warn_str = n.name * " on tree but not found in names."
            if tolerate_missing == 0
                @error warn_str
            end
            if tolerate_missing == 1
                @warn warn_str
            end
            uninformative_message!(n.message)
        end
    end
end

function internal_message_init!(tree::FelNode, empty_message::LazyPartition{T}) where {T <: Partition}
    for node in getnodelist(tree)
        if !isleafnode(node)
            node.child_messages = []
            for _ in node.children
                child_message = LazyPartition{T}(nothing)
                child_message.memoryblocks = empty_message.memoryblocks
                push!(node.child_messages, [child_message])
            end
        end

        message = LazyPartition{T}(nothing)
        parent_message = LazyPartition{T}(nothing)
        message.memoryblocks = empty_message.memoryblocks
        parent_message.memoryblocks = empty_message.memoryblocks
        node.message = [message]
        node.parent_message = [parent_message]
    end
end

#Need a faster realization of the leafnodes that avoids allocating new matrices and vectors...
function lazy_obs2partition!(dest::LazyPartition{PType}, seq) where {PType}
end