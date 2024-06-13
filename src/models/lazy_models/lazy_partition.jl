export LazyPartition
"""
# Constructor
    LazyPartition{PType}()
Initialize an empty `LazyPartition` that is meant for wrapping a partition of type `PType`.

# Description
With this data structure, you can wrap a partition of choice. 
The idea is that in some message passing algorithms, there is only a wave of partitions which need to actualize. 
For instance, a wave following a root-leaf path, or a depth-first traversal.
In which case, we can be more economical with our memory consumption.
With a worst case memory complexity of O(log(n)), where n is the number of nodes, functionality is provided for:
- `log_likelihood!`
- `felsenstein!`
- `sample_down!`
!!! note
    For successive `felsenstein!` calls, we need to extract the information at the root somehow after each call. This can be done with e.g. `total_LL` or `site_LLs`.

# Further requirements
Suppose you want to wrap a partition of `PType` with `LazyPartition`:
- If you're calling `log_likelihood!` and `felsenstein!`:
    - `obs2partition!(partition::PType, obs)` that transforms an observation to a partition.
- If you're calling `sample_down!`:
    - `partition2obs(partition::PType)` that returns the most likely state from a partition, inverts `obs2partition!`.
"""
mutable struct LazyPartition{PType} <: Partition where {PType <: Partition}
    partition::Union{PType, Nothing}
    memoryblocks::Vector{PType} #A stack that is meant for LazyPartitions in the same tree to share
    static::Bool #Does the LazyPartition keep its partition during a message passing wave? Defaults to false
    obs #Leaf nodes can store their observations here in a preferably compact data structure. Note: allowing Any type is probably bad for performance

    function LazyPartition{PType}(partition, memoryblocks) where {PType <: Partition}
        new(partition, memoryblocks, false)
    end
end

function LazyPartition{PType}() where {PType <: Partition}
    LazyPartition{PType}(nothing, Vector{PType}())
end

export LazyUp
"""
# Constructor
    LazyUp()

# Description
Indicate that we want to do an upward pass, e.g. `felsenstein!`.
"""
struct LazyUp <: LazyDirection end
export LazyDown
"""
# Constructors
    LazyDown(stores_obs)
    LazyDown() = LazyDown(x::FelNode -> true)

# Description
Indicate that we want to do a downward pass, e.g. `sample_down!`.
The function passed to the constructor takes a `node::FelNode` as input and returns a `Bool` that decides if `node` stores its observations.
"""
struct LazyDown <: LazyDirection
    stores_obs #Function: FelNode -> Bool. Does the node store its observation after a down pass?
    
    function LazyDown(stores_obs)
        new(stores_obs)
    end
end

LazyDown() = LazyDown(x::FelNode -> true) #Every node stores obs

"""
I store:
    - References to partitions in a global stack
    - Observations in the obs field of the LazyPartition

To achieve a maximum amount of memoryblocks available during felsenstein!, the workflow of backward! and combine! is this:
    - We pop a memoryblock of the stack and put it in dest.partition prior to a backward!.
    Now, we perform the backward!, then push our source.partition to the memoryblock stack
    - After a combine!, we push src.partition to the memoryblock stack
"""

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
    (isnothing(src.partition) || isnothing(dest.partition)) && throw(ArgumentError("The partition field in both the source and destination LazyPartitions must be defined to combine them.\nNote that LazyPartition can only be used for an upward Felsenstein pass or sample_down!."))
    combine!(dest.partition, src.partition)
    safe_release_partition!(src)
end


"""
The idea with backward!{LazyPartition} is to avoid having to call the GC too often,
therefore, we need access to shared memoryblocks during the message passing wave. This is only really necessary if PType <: MultisitePartition,
but for consistency, we might as well do it for all types?
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
        obs2partition!(source.partition, source.obs)
    end
    dest.partition = pop!(source.memoryblocks)
    backward!(dest.partition, source.partition, model, node)
    safe_release_partition!(source)
end

function forward!(
    dest::LazyPartition{PType},
    source::LazyPartition{PType},
    model::BranchModel,
    node::FelNode,
) where {PType <: Partition}
    isnothing(source.partition) && throw(ArgumentError("The partition field in the source LazyPartition must be defined for a forward! call.\nNote that LazyPartition can only be used for an upward Felsenstein pass or sample_down!."))
    if isnothing(dest.partition)
        isempty(dest.memoryblocks) && throw(ArgumentError("There needs to be atleast 1 available memoryblock for a forward!.\nNote that LazyPartition can only be used for an upward Felsenstein pass or sample_down!."))
        dest.partition = pop!(dest.memoryblocks)
    end
    if !isroot(node) && node == node.parent.children[end] #We no longer need source to be static
        source.static = false #Maybe this is something we only want to do during sample_down!
    end
    forward!(dest.partition, source.partition, model, node)
    safe_release_partition!(source)
end

#General note: after a felsenstein!, this must be called on root.message for it to be recycled into memoryblocks.
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
- Should be run on a tree containing LazyPartitions before running `felsenstein!`. Sorts for a minimal count of active partitions during a felsenstein!
- Returns the minimum length of memoryblocks (-1) required for a `felsenstein!` prop. We need a temporary memoryblock during `backward!`, hence the '-1'.
!!! note
    Since felsenstein! uses a stack, we want to avoid having long node.children[1].children[1]... chains
"""
function lazysort!(tree::FelNode)
    node_stack = [(tree, true)]
    result_stack = Vector{Int64}() #Result is the maximum active partitions a node requires
    while !isempty(node_stack)
        node, first = pop!(node_stack)
        if !isleafnode(node)
            if first
                push!(node_stack, (node, false))
                for child in node.children
                    push!(node_stack, (child, true))
                end
            end
            if !first
                #maximum_active_partitions maps directly onto node.children since we have two pop! operations (reverses order) that cancel out
                maximum_active_partitions = [pop!(result_stack) for _ in node.children]
                perm = sortperm(maximum_active_partitions)
                node.children = node.children[perm] #Sort to avoid long node.children[1].children[1]... chains
                for (i, idx) in enumerate(perm)
                    maximum_active_partitions[idx] += length(node.children) - i #Add n.o. nodes that is to the right
                end
                push!(result_stack, maximum(maximum_active_partitions))
            end
        else
            push!(result_stack, 1)
        end
    end
    return pop!(result_stack)
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
    lazyprep!(tree::FelNode, initial_message::Vector{<:Partition}; partition_list = 1:length(tree.message), direction::LazyDirection = LazyUp())

Extra, intermediate step of tree preparations between initializing messages across the tree and calling message passing algorithms with `LazyPartition`.
1. Perform a `lazysort!` on tree to obtain the optimal tree for a lazy `felsenstein!` prop, or a `sample_down!`.
2. Fix `tree.parent_message` to an initial message.
3. Preallocate sufficiently many inner partitions needed for a `felsenstein!` prop, or a `sample_down!`.
4. Specialized preparations based on the direction of the operations (`forward!`, `backward!`). `LazyDown` or `LazyUp`.

See also `LazyDown`, `LazyUp`.
"""
function lazyprep!(tree::FelNode, initial_message::Vector{<:Partition}; partition_list = 1:length(tree.message), direction::LazyDirection = LazyUp())
    @assert length(initial_message) == length(partition_list)
    maximum_active_partitions = lazysort!(tree) + 1 # the +1 comes from using an extra temporary memoryblock during backward!
    for (p, part) in enumerate(partition_list)
        parent_part = tree.parent_message[part]
        parent_part.partition = initial_message[p]
        parent_part.static = true
        for _ = 1:maximum_active_partitions
            push!(parent_part.memoryblocks, partition_from_template(initial_message[p]))
        end
    end
    lazyprep!(tree, direction, partition_list = partition_list) #Specialized prep that dispatches on direction
    return maximum_active_partitions
end

lazyprep!(
    tree::FelNode, 
    initial_partition::Partition; 
    partition_list = 1:length(tree.message), 
    direction::LazyDirection = LazyUp()
) = lazyprep!(tree, [initial_partition], partition_list = partition_list, direction = direction)


#The reason we can have the same sort for felsenstein! and sample_down! is that two negatives cancels out. Upward vs. Downward + Recursive vs. Iterative (stack)
function sample_partition!(part::LazyPartition)
    sample_partition!(part.partition)
    if !isdefined(part, :obs)
        return
    end
    part.obs = partition2obs(part.partition)
    safe_release_partition!(part)
end

function lazyprep!(tree::FelNode, direction::LazyUp; partition_list = 1:length(tree.message))
    return
end

function lazyprep!(tree::FelNode, direction::LazyDown; partition_list = 1:length(tree.message))
    for node in getnodelist(tree)
        if direction.stores_obs(node)
            property = :obs
            value = nothing #Signal that something should end up here
            setproperty!.(node.message[partition_list], property, value)
        end
        if !isleafnode(node)
            property = :static
            value = true #Static until finished with all forwarding as source
            setproperty!.(node.message[partition_list], property, value)
        end
    end
end
