struct InternalStandardRootOpt <: RootOpt
    root_LL!::Function
    K::Int
end

Base.length(root_opt::InternalStandardRootOpt) = root_opt.K

mutable struct StandardRootOpt <: RootOpt
    K::Int
end

mutable struct StandardRootSample <: UniformRootPositionSample
    acc_ratio::Tuple{Float64, Int64, Int64}
    radius::Float64
    consecutive::Int64

    function StandardRootSample(radius::Float64,consecutive::Int64)
        new((0.0, 0, 0), radius, consecutive)
    end
end

Base.length(root_sample::StandardRootSample) = root_sample.consecutive
radius(root_sample::StandardRootSample) = root_sample.radius

#Assume that felsenstein_roundtrip! has been called
#Compute the log likelihood of observations below this root-candidate
function root_LL_below!(
    dest::Vector{<:Partition},
    temp::Vector{<:Partition},
    dist_above_node::Real,
    node::FelNode,
    model_list::Vector{<:BranchModel};
    partition_list = 1:length(tree.message)
)
    @assert 0.0 <= dist_above_node < node.branchlength || dist_above_node == node.branchlength == 0.0 #if dist_above_node == node.branchlength != 0.0, then it's node.parent with 0.0 dist_above_child that should be called
    branchlength = node.branchlength
    for (p, part) in enumerate(partition_list)
        node.branchlength = dist_above_node
        backward!(dest[p], node.message[part], model_list[part], node)
        node.branchlength = branchlength - dist_above_node
        #We don't want to double count parent message of root
        if !(isroot(node) && dist_above_node == zero(dist_above_node))
            backward!(temp[p], node.parent_message[part], model_list[part], node)
            combine!(dest[p], temp[p])
        end
    end
    node.branchlength = branchlength
end

function steal_messages!(new_root::FelNode, old_root::FelNode)
    new_root.message = old_root.message
    new_root.parent_message = old_root.parent_message
    new_root.child_messages = old_root.child_messages
end

function default_root_LL_wrapper(parent_message::Vector{<:Partition})
    function root_LL!(message::Vector{<:Partition})
        combine!.(message, parent_message)
        return parent_message, sum(total_LL.(message))
    end
    return root_LL!
end

root_update!(
    root_update::RootUpdate,
    tree::FelNode,
    models::Vector{<:BranchModel};
    kwargs...
) = root_update!(
        root_update,
        tree,
        x -> models;
        kwargs...
    )

root_update!(
    root_update::RootUpdate,
    tree::FelNode,
    model::BranchModel;
    kwargs...
) = root_update!(
        root_update,
        tree,
        x -> [model];
        kwargs...
    )

"""
    root_update!(root_update::RootUpdate, tree::FelNode, models; partition_list = 1:length(tree.message))

A more general version of [`root_optim!`](@ref). Here `root_update` can be either an optimization or a sampling (or more generally, a RootUpdate).
"""
function root_update!(root_update::RootUpdate, tree::FelNode, models; partition_list = 1:length(tree.message))
    #Initialize some messages
    node_message = copy_message(tree.parent_message[partition_list])
    temp_message = copy_message(tree.parent_message[partition_list])

    #Optimize the root position + root state
    new_value = root_update(tree, models, partition_list, node_message, temp_message)

    #Update the root position + root state
    new_root = new_value.root == tree ? tree : reroot!(new_value.root, dist_above_child = new_value.dist_above_node) #Maybe reroot! should take care of this?
    steal_messages!(new_root, tree)
    new_root.parent_message[partition_list] .= new_value.state
    return new_root
end

"""
    root_optim!(tree::FelNode, models; <keyword arguments>)

Optimizes the root position and root state of a tree. Returns the new, optimal root node.
models can either be a single model (if the messages on the tree contain just one Partition) or an array of models, if the messages have >1 Partition, or 
a function that takes a node, and returns a Vector{<:BranchModel} if you need the models to vary from one branch to another.

# Keyword Arguments
- `partition_list=1:length(tree.message)`: (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over (but you probably want to optimize root position and root state with all models, the default option).
- `root_LL!=default_root_LL_wrapper(tree.parent_message[partition_list])`: a function that takes a message and returns a (optimal) parent message and LL (log likelihood). The default option uses the constant `tree.parent_message[partition_list]` as parent message for all root-candidates.
- `K=10`: the number of equidistant root-candidate points along a branch. (only to be used in the frequentist framework!?)
"""
function root_optim!(tree::FelNode, models; partition_list = 1:length(tree.message), root_LL! = default_root_LL_wrapper(tree.parent_message[partition_list]), K = 10)
    return root_update!(InternalStandardRootOpt(root_LL!, K), tree, models, partition_list = partition_list)
end

function (root_opt::RootOpt)(tree::FelNode, models, partition_list, node_message::Vector{<:Partition}, temp_message::Vector{<:Partition})

    #Do most of the message passing
    felsenstein_roundtrip!(tree, models, partition_list = partition_list, temp_message = temp_message)

    #Initialize the fallback optimum
    opt_root = tree
    opt_dist = 0.0
    opt_LL = log_likelihood(tree, models, partition_list = partition_list)
    opt_starting_message = copy_message(tree.parent_message[partition_list])

    #Optimize the root position + root state
    nodelist = getnodelist(tree)
    for node in nodelist
        model_list = models(node)
        for dist_above_node in unique(range(0.0, node.branchlength, length=length(root_opt)+1)[1:end-1])
            #                  unique() to avoid recomputations
            #Compute the log likelihood of observations below this root-candidate...
            root_LL_below!(
                node_message,
                temp_message,
                dist_above_node,
                node,
                model_list,
                partition_list = partition_list
            )
            node_starting_message, LL = root_opt(node_message)
            if LL > opt_LL
                opt_root, opt_dist, opt_LL = node, dist_above_node, LL
                copy_partition_to!.(opt_starting_message, node_starting_message)
            end
        end
    end
    return (root=opt_root, dist_above_node=opt_dist, state=opt_starting_message)
end

#(opt::RootOpt)(tree::FelNode, model::BranchModel, args...; kwargs...) = opt(tree, x -> [model], args...; kwargs...)

#(opt::RootOpt)(tree::FelNode, models::Vector{<:BranchModel}, args...; kwargs...) = opt(tree, x -> models, args...; kwargs...)




(root_opt::InternalStandardRootOpt)(message::Vector{<:Partition}) = root_opt.root_LL!(message)

function (root_opt::StandardRootOpt)(
    tree::FelNode,
    models,
    partition_list,
    node_message::Vector{<:Partition},
    temp_message::Vector{<:Partition},
)
    InternalStandardRootOpt(
        default_root_LL_wrapper(tree.parent_message[partition_list]),
        root_opt.K,
    )(
        tree,
        models,
        partition_list,
        node_message,
        temp_message,
    )
end

(root_opt::StandardRootOpt)(tree::FelNode, model::BranchModel, partition_list, node_message::Vector{<:Partition}, temp_message::Vector{<:Partition}) = root_opt(tree, x -> [model], partition_list, node_message, temp_message)


(root_opt::StandardRootOpt)(tree::FelNode, models::Vector{<:BranchModel}, partition_list, node_message::Vector{<:Partition}, temp_message::Vector{<:Partition}) = root_opt(tree, x -> models, partition_list, node_message, temp_message)



function (root_sample::RootSample)(tree::FelNode, models, partition_list, node_message::Vector{<:Partition}, temp_message::Vector{<:Partition})
    #Do most of the message passing
    felsenstein_roundtrip!(tree, models, partition_list = partition_list, temp_message = temp_message)

    #Sample the root position + root state
    function LL!(position::@NamedTuple{root::FelNode, dist_above_node::Float64}, state::Vector{<:Partition})
        root_LL_below!(node_message, temp_message, position.dist_above_node, position.root, models(position.root), partition_list = partition_list)
        combine!.(node_message, state)
        return sum(total_LL.(node_message))
    end
    #Initialize variables
    sampled_position = (root=tree, dist_above_node=0.0)
    sampled_state = tree.parent_message[partition_list]
    for _ = 1:length(root_sample) #this solution is arguably hacky, but we want to squeeze out as much of `felsenstein_roundtrip!` as possible
        sampled_position = metropolis_step(Base.Fix2(LL!, sampled_state), root_sample, sampled_position)
        sampled_state = metropolis_step(Base.Fix1(LL!, sampled_position), root_sample, sampled_state)
    end
    
    return merge(sampled_position, (state=sampled_state,))
end


function traverse(node0::FelNode, dist_above_node0::Float64, radius::Float64)
    stack = [(node0, dist_above_node0, radius, length(node0.children)+1)]
    points = Vector{Tuple{FelNode, Float64, Int64}}()
    weights = Vector{Float64}()

    while !isempty(stack)
        node, dist_above_node, radius_left, prev_ind = pop!(stack)
                                            #interpret 1,2,... as which child ind we're coming from (if it's out of range, we're at node0), 
                                            #0 as we're coming from root,
        if prev_ind > 0 && !isroot(node) #Upward traversal
            child_ind = findfirst(x -> x == node, node.parent.children)
            radius_that_would_be_left = radius_left - (node.branchlength - dist_above_node)
            if radius_that_would_be_left < 0.0
                push!(points, (node, dist_above_node + radius_left, child_ind))
                push!(weights, radius_left)
            else
                radius_that_would_be_left != 0.0 && push!(stack, (node.parent, 0.0, radius_that_would_be_left, child_ind))
                push!(points, (node, node.branchlength, child_ind)) #(node, node.branchlength) <=> (node.parent, 0.0)
                push!(weights, node.branchlength - dist_above_node)
            end
        end
        #Downward traversal
        radius_that_would_be_left = radius_left - dist_above_node
        if radius_that_would_be_left < 0.0
            push!(points, (node, -radius_that_would_be_left, 0))
            push!(weights, radius_left)
            continue
        end
        if dist_above_node != 0.0 #node0 may add itself, otherwise dist_above_node is 0.0
            push!(points, (node, 0.0, 0)) #most cases down will be true, but when it is false, we want to remember that
            push!(weights, dist_above_node)
        end
        for (child_ind, child) in enumerate(node.children)
            if child_ind == prev_ind
                continue
            end
            radius_that_would_be_left = radius_left - dist_above_node - child.branchlength
            if radius_that_would_be_left < 0.0
                push!(points, (child, -radius_that_would_be_left, 0))
                push!(weights, radius_left - dist_above_node)
                continue
            end
            radius_that_would_be_left != 0.0 && push!(stack, (child, 0.0, radius_that_would_be_left, 0))
            push!(points, (child, 0.0, 0))
            push!(weights, child.branchlength)
        end
    end
    return points, weights
end

function log_proposal(modifier::UniformRootPositionSample,
    x::@NamedTuple{root::FelNode, dist_above_node::Float64},
    conditioned_on::@NamedTuple{root::FelNode, dist_above_node::Float64})
    points, weights = traverse(conditioned_on..., radius(modifier))
    #if x is not within radius(modifier) of conditioned_on, then we should return log(0.0)
    return -log(sum(weights))
end

# Propose a new root position with a local uniform distribution
function proposal(modifier::UniformRootPositionSample, curr_value::@NamedTuple{root::FelNode, dist_above_node::Float64})
    points, weights = traverse(curr_value..., radius(modifier))
    #Sample a new root position within radius radius(modifier)
    cum = cumsum(weights)
    sample = rand() * cum[end]
    idx = searchsortedfirst(cum, sample)
    diff = cum[idx] - sample
    node, dist_above_node, prev_ind = points[idx]
    dist_above_node += prev_ind > 0 ? -diff : diff
    return (root=node, dist_above_node=dist_above_node)
end

log_prior(::UniformRootPositionSample, curr_value::@NamedTuple{root::FelNode, dist_above_node::Float64}) = 0.0 #Uninformative/improper prior

function proposal(::StandardRootSample, curr_value::Vector{<:Partition})
    return curr_value
end

log_prior(::StandardRootSample, curr_value::Vector{<:Partition}) = 0.0 #Uninformative/improper prior

#Note: I've got a feeling there's no difference in expressibility between RootOpt and StandardRootOpt, if (::<:RootOpt)(tree::FelNode, ...) is not overloaded. 