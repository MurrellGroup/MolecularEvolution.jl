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

"""
    root_optim!(tree::FelNode, models; <keyword arguments>)

Optimizes the root position and root state of a tree. Returns the new, optimal root node.
models can either be a single model (if the messages on the tree contain just one Partition) or an array of models, if the messages have >1 Partition, or 
a function that takes a node, and returns a Vector{<:BranchModel} if you need the models to vary from one branch to another.

# Keyword Arguments
- `partition_list=1:length(tree.message)`: (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over (but you probably want to optimize root position and root state with all models, the default option).
- `starting_message_modifier!=(objective, starting_message0::Vector{<:Partition}) -> starting_message0`: (can either be) an optimizer (or a sampler) of the root state. `objective` returns the log likelihood of a root state (and implicitly, a root position).
- `starting_message0::Vector{<:Partition}=copy_message(tree.parent_message[partition_list])`: the initial starting message used by `starting_message_modifier!`.
- `K=10`: the number of equidistant root-candidate points along a branch. (only to be used in the frequentist framework!?)
"""
function root_optim!(
    tree::FelNode,
    models;
    partition_list = 1:length(tree.message),
    starting_message_modifier! = (objective, starting_message0::Vector{<:Partition}) -> starting_message0,
    starting_message0::Vector{<:Partition} = copy_message(tree.parent_message[partition_list]),
    K = 10 #Number of root-candidate points on a branch
)
    #Initialize some messages
    node_message = copy_message(tree.parent_message[partition_list])
    node_starting_message = copy_message(tree.parent_message[partition_list])
    temp_message = copy_message(tree.parent_message[partition_list])

    #Initialize the fallback optimum
    opt_root = tree
    opt_dist = 0.0
    opt_LL = log_likelihood(tree, models, partition_list = partition_list)
    opt_starting_message = copy_message(tree.parent_message[partition_list])

    #Do most of the message passing
    felsenstein_roundtrip!(tree, models, partition_list = partition_list, temp_message = temp_message)

    #Optimize the root position + root state
    nodelist = getnodelist(tree)
    for node in nodelist
        copy_partition_to!.(node_starting_message, starting_message0)
        model_list = models(node)
        for dist_above_node in unique(range(0.0, node.branchlength, K + 1)[1:end-1])
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
            function objective(starting_message::Vector{<:Partition})
                for p = 1:length(partition_list)
                    #... combine it with a root state...
                    combine!(temp_message[p], [starting_message[p], node_message[p]], true)
                end
                #... and get the total log likelihood.
                return sum(total_LL.(temp_message))
            end
            node_starting_message = starting_message_modifier!(objective, node_starting_message)
            #Reuse this as the starting_message0 for the next iteration
            LL = objective(node_starting_message)
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
    starting_message_modifier! = (objective, starting_message0::Vector{<:Partition}) -> starting_message0,
    starting_message0::Vector{<:Partition} = copy_message(tree.parent_message[partition_list]),
    K = 10 #Number of root-candidate points on a branch
) = root_optim!(
        tree,
        x -> models,
        partition_list = partition_list,
        starting_message_modifier! = starting_message_modifier!,
        starting_message0 = starting_message0,
        K = K
    )

root_optim!(
    tree::FelNode,
    model::BranchModel;
    partition_list = 1:length(tree.message),
    starting_message_modifier! = (objective, starting_message0::Vector{<:Partition}) -> starting_message0,
    starting_message0::Vector{<:Partition} = copy_message(tree.parent_message[partition_list]),
    K = 10 #Number of root-candidate points on a branch
) = root_optim!(
        tree,
        x -> [model],
        partition_list = partition_list,
        starting_message_modifier! = starting_message_modifier!,
        starting_message0 = starting_message0,
        K = K
    )