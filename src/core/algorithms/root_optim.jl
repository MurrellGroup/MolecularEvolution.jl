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
        backward!(temp[p], node.parent_message[part], model_list[part], node)
        combine!(dest[p], temp[p])
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
function root_optim!(
    tree::FelNode,
    models;
    partition_list = 1:length(tree.message),
    root_LL! = default_root_LL_wrapper(tree.parent_message[partition_list]),
    K = 10 #Number of root-candidate points on a branch
)
    #Initialize some messages
    node_message = copy_message(tree.parent_message[partition_list])
    temp_message = copy_message(tree.parent_message[partition_list])

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
        for dist_above_node in unique(range(0.0, node.branchlength, length=K+1)[1:end-1])
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
            node_starting_message, LL = root_LL!(node_message)
            #TODO: enable root sampling
            if LL > opt_LL
                opt_root, opt_dist, opt_LL = node, dist_above_node, LL
                copy_partition_to!.(opt_starting_message, node_starting_message)
            end
        end
    end
    new_root = opt_root == tree ? tree : reroot!(opt_root, dist_above_child = opt_dist) #Maybe reroot! should take care of this?
    steal_messages!(new_root, tree)
    new_root.parent_message[partition_list] .= opt_starting_message
    return new_root
end

root_optim!(
    tree::FelNode,
    models::Vector{<:BranchModel};
    partition_list = 1:length(tree.message),
    root_LL! = default_root_LL_wrapper(tree.parent_message[partition_list]),
    K = 10 #Number of root-candidate points on a branch
) = root_optim!(
        tree,
        x -> models,
        partition_list = partition_list,
        root_LL! = root_LL!,
        K = K
    )

root_optim!(
    tree::FelNode,
    model::BranchModel;
    partition_list = 1:length(tree.message),
    root_LL! = default_root_LL_wrapper(tree.parent_message[partition_list]),
    K = 10 #Number of root-candidate points on a branch
) = root_optim!(
        tree,
        x -> [model],
        partition_list = partition_list,
        root_LL! = root_LL!,
        K = K
    )

function root_position_sample!(
    tree::FelNode,
    models;
    partition_list = 1:length(tree.message),
    root_LL! = default_root_LL_wrapper(tree.parent_message[partition_list]),
    root_sampler::UnivariateSampler = RootPositionSampler(),
    K = 1 #Number of consecutive metropolis steps
)
    #Initialize some messages
    node_message = copy_message(tree.parent_message[partition_list])
    temp_message = copy_message(tree.parent_message[partition_list])

    #Do most of the message passing
    felsenstein_roundtrip!(tree, models, partition_list = partition_list, temp_message = temp_message)

    #Initialize variables
    sampled_root = tree
    sampled_dist = 0.0
    starting_message = copy_message(tree.parent_message[partition_list])

    #Sample the root position + root state
    function LL_position!(curr_value::@NamedTuple{root::FelNode, dist_above_node::Real})
        root_LL_below!(node_message, temp_message, curr_value.dist_above_node, curr_value.root, models(curr_value.root), partition_list = partition_list)
        combine!.(node_message, temp_message)
        return sum(total_LL.(node_message))
    end
    new_value = metropolis_step(LL_position!, root_sampler, (root=tree, dist_above_node=0.0))
    sampled_root, sampled_dist = new_value.root, new_value.dist_above_node
    new_root = sampled_root == tree ? tree : reroot!(sampled_root, dist_above_child = sampled_dist) #Maybe reroot! should take care of this?
    steal_messages!(new_root, tree)
    new_root.parent_message[partition_list] .= starting_message
    return new_root
end

struct RootPositionSampler <: UnivariateSampler
    acc_ratio::Vector{Int}
    function RootPositionSampler()
        new([0,0])
    end
end 

# Propose a new root position with a global uniform distribution
function proposal(sampler::RootPositionSampler, curr_value::@NamedTuple{root::FelNode, dist_above_node::Real})
    nodelist = getnodelist(curr_value.root)
    cum = cumsum(n.branchlength for n in nodelist)
    sample = rand() * cum[end]
    idx = searchsortedfirst(cum, sample)
    return (nodelist[idx], cum[idx] - sample)
end

log_prior(sampler::RootPositionSampler, curr_value) = 0.0 #Uninformative/improper prior


#=TODO: fix this mess
-We want to conform to the root_update!(dispatch_thing, tree, models) pattern
-I need to figure out a generic root_sample! and not only root_position_sample!, in which case, I have to think about if I want to call 

---------
#There isn't a straightforward generic root_update! yet, therefore no abstract type RootModifier
abstract type RootOpt end
abstract type RootSampler end
=#
abstract type RootUpdate <: Function end
abstract type RootOpt <: RootUpdate end
abstract type RootSample <: RootUpdate end

struct StandardRootOpt <: RootOpt
    root_LL!::Function
end

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

function root_optim!(tree::FelNode, models; root_LL! = default_root_LL_wrapper(tree.parent_message[partition_list]), kwargs...)
    return root_update!(StandardRootOpt(root_LL!), tree, models; kwargs...)
end

function (root_opt::RootOpt)(tree::FelNode, models, partition_list, node_message::Vector{<:Partition}, temp_message::Vector{<:Partition}; K = 10)

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
        for dist_above_node in unique(range(0.0, node.branchlength, length=K+1)[1:end-1])
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

(root_opt::StandardRootOpt)(message::Vector{<:Partition}) = root_opt.root_LL!(message)

function (root_sample::RootSample)(tree::FelNode, models, partition_list, node_message::Vector{<:Partition}, temp_message::Vector{<:Partition})
    #Do most of the message passing
    felsenstein_roundtrip!(tree, models, partition_list = partition_list, temp_message = temp_message)

    #Initialize variable
    starting_message = copy_message(tree.parent_message[partition_list])

    #Sample the root position + root state
    function LL!(position::@NamedTuple{root::FelNode, dist_above_node::Real}, state::Vector{<:Partition})
        root_LL_below!(node_message, temp_message, position.dist_above_node, position.root, models(position.root), partition_list = partition_list)
        combine!.(node_message, state)
        return sum(total_LL.(node_message))
    end
    sampled_position = metropolis_step(Base.Fix2(LL!, starting_message), root_sampler, (root=tree, dist_above_node=0.0))
    sampled_state = metropolis_step(Base.Fix1(LL!, sampled_position), root_sampler, starting_message)
    return merge(sampled_position, (state=sampled_state,))
end

struct UniformRootPositionSampler <: RootSample
    acc_ratio::Vector{Int}
    function UniformRootPositionSampler()
        new([0,0])
    end
end 

# Propose a new root position with a global uniform distribution
function proposal(sampler::UniformRootPositionSampler, curr_value::@NamedTuple{root::FelNode, dist_above_node::Real})
    nodelist = getnodelist(curr_value.root)
    cum = cumsum(n.branchlength for n in nodelist)
    sample = rand() * cum[end]
    idx = searchsortedfirst(cum, sample)
    return (nodelist[idx], cum[idx] - sample)
end

function proposal(sampler::UniformRootPositionSampler, curr_value::Vector(<:Partition))
    return curr_value
end

log_prior(sampler::RootPositionSampler, curr_value) = 0.0 #Uninformative/improper prior