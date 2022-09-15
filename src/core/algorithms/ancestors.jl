#BM: Merged a number of different things into the same framework.
#There are a family of things we want to do over ancestral states that all have the following form:
#1. Felsenstein up pass
#2. Recurse down the tree, doing:
#3. node reconstruction = R(node, parent reconstruction)
#   where R is a function that takes a message and returns a "reconstructed" message.
#We want to do this without disturbing the up and down messages, so that you get the same answer when you re-run.

#This might have a better home in the general tree utils section.
export depth_first_traversal
function depth_first_traversal(root, func)
    stack = [root]
    while length(stack) > 0
        node = pop!(stack)
        func(node)
        for child in node.children
            push!(stack, child)
        end
    end
end


export depth_first_reconstruction
#I think this is the right one to overload, with the model_func vs model vs models
function depth_first_reconstruction(
    tree,
    r,
    model_func;
    run_fel_up = true,
    run_fel_down = true,
    partition_list = 1:length(tree.message),
    node_message_dict = Dict{FelNode,Vector{Partition}}(),
)
    if run_fel_up
        felsenstein!(tree, model_func, partition_list = partition_list)
    end
    if run_fel_down
        felsenstein_down!(tree, model_func, partition_list = partition_list)
    end
    stack = [tree]
    while length(stack) > 0
        node = pop!(stack)
        models = model_func(node)
        r(node_message_dict, node, models, partition_list)
        for child in node.children
            push!(stack, child)
        end
    end
    return node_message_dict
end

function depth_first_reconstruction(
    tree,
    r,
    model::BranchModel;
    run_fel_up = true,
    run_fel_down = true,
    partition_list = 1:length(tree.message),
    node_message_dict = Dict{FelNode,Vector{Partition}}(),
)
    depth_first_reconstruction(
        tree,
        r,
        x -> [model],
        run_fel_up = run_fel_up,
        run_fel_down = run_fel_down,
        partition_list = partition_list,
        node_message_dict = node_message_dict,
    )
end

function depth_first_reconstruction(
    tree,
    r,
    models::Vector{<:BranchModel};
    run_fel_up = true,
    run_fel_down = true,
    partition_list = 1:length(tree.message),
    node_message_dict = Dict{FelNode,Vector{Partition}}(),
)
    depth_first_reconstruction(
        tree,
        r,
        x -> models,
        run_fel_up = run_fel_up,
        run_fel_down = run_fel_down,
        partition_list = partition_list,
        node_message_dict = node_message_dict,
    )
end

#For marginal reconstructions
function reconstruct_marginal_node!(
    node_message_dict::Dict{FelNode,Vector{Partition}},
    node::FelNode,
    model_array::Vector{<:BranchModel},
    partition_list,
)
    #Note: these family of functions can be passed a dictionary for repeat use, but then the "structure" must be unchanged.
    if !haskey(node_message_dict, node)
        node_message_dict[node] = deepcopy(node.message[partition_list])
    end
    m = node_message_dict[node]
    for (p, part) in enumerate(partition_list)
        forward!(m[p], node.parent_message[part], model_array[part], node)
        combine!(m[p], node.message[part])
    end
end

export marginal_state_dict
"""
    marginal_state_dict(tree::FelNode, model; partition_list = 1:length(tree.message), node_message_dict = Dict{FelNode,Vector{Partition}}())

Takes in a tree and a model (which can be a single model, an array of models, or a function that maps FelNode->Array{<:BranchModel}), and
returns a dictionary mapping nodes to their marginal reconstructions (ie. P(state|all observations,model)). A subset of partitions can be specified by partition_list,
and a dictionary can be passed in to avoid re-allocating memory, in case you're running this over and over.
"""
function marginal_state_dict(
    tree::FelNode,
    model;
    partition_list = 1:length(tree.message),
    node_message_dict = Dict{FelNode,Vector{Partition}}(),
)
    return depth_first_reconstruction(
        tree,
        reconstruct_marginal_node!,
        model,
        partition_list = partition_list,
        node_message_dict = node_message_dict,
    )
end

#For joint max reconstructions
export dependent_reconstruction!
function dependent_reconstruction!(
    node_message_dict::Dict{FelNode,Vector{Partition}},
    node::FelNode,
    model_array::Vector{<:BranchModel},
    partition_list;
    f! = x -> nothing,
)
    if !haskey(node_message_dict, node)
        node_message_dict[node] = deepcopy(node.message[partition_list])
    end
    m = node_message_dict[node]
    for (p, part) in enumerate(partition_list)
        if isroot(node)
            forward!(m[p], node.parent_message[part], model_array[part], node)
        else
            forward!(m[p], node_message_dict[node.parent][p], model_array[part], node)
        end
        combine!(m[p], node.message[part])
        f!(m[p])
    end
end


#BM: This algorithm starts from the standard felsenstein upward pass, and then selects the most likely state at the root.
#Then, conditioned on that state, and the upward pass of the children of the root, it selects the max of each child.
#This is something like: the maximum likelihood of internal node states, conditioned on the maximum marginal likelihood of the root state.
#An interesting case to think about is when you're rooted on a leaf node, so the max marginal and max joint at the root are identical.
#Then probs from that node will always incur a "max_partition!" call after combining with the upward pass to that node.
#So where does that leave us?

#For joint max reconstructions
reconstruct_cascading_max_node!(node_message_dict, node, model_array, partition_list) =
    dependent_reconstruction!(
        node_message_dict,
        node,
        model_array,
        partition_list,
        f! = max_partition!,
    )
export cascading_max_state_dict
"""
    cascading_max_state_dict(tree::FelNode, model; partition_list = 1:length(tree.message), node_message_dict = Dict{FelNode,Vector{Partition}}())

Takes in a tree and a model (which can be a single model, an array of models, or a function that maps FelNode->Array{<:BranchModel}), and
returns a dictionary mapping nodes to their inferred ancestors under the following scheme: the state that maximizes the marginal likelihood is selected at the root,
and then, for each node, the maximum likelihood state is selected conditioned on the maximized state of the parent node and the observations of all descendents.
A subset of partitions can be specified by partition_list, and a dictionary can be passed in to avoid re-allocating memory, in case you're running this over and over.
"""
function cascading_max_state_dict(
    tree::FelNode,
    model;
    partition_list = 1:length(tree.message),
    node_message_dict = Dict{FelNode,Vector{Partition}}(),
)
    return depth_first_reconstruction(
        tree,
        reconstruct_cascading_max_node!,
        model,
        partition_list = partition_list,
        node_message_dict = node_message_dict,
    )
end

#For endpoint conditioned sampling - could consider merging this into the above function
conditioned_sample_node!(node_message_dict, node, model_array, partition_list) =
    dependent_reconstruction!(
        node_message_dict,
        node,
        model_array,
        partition_list,
        f! = sample_partition!,
    )
export endpoint_conditioned_sample_state_dict
"""
    endpoint_conditioned_sample_state_dict(tree::FelNode, model; partition_list = 1:length(tree.message), node_message_dict = Dict{FelNode,Vector{Partition}}())

Takes in a tree and a model (which can be a single model, an array of models, or a function that maps FelNode->Array{<:BranchModel}), and draws samples under the model
conditions on the leaf observations. These samples are stored in the node_message_dict, which is returned. A subset of partitions can be specified by partition_list, and a
dictionary can be passed in to avoid re-allocating memory, in case you're running this over and over.
"""
function endpoint_conditioned_sample_state_dict(
    tree::FelNode,
    model;
    partition_list = 1:length(tree.message),
    node_message_dict = Dict{FelNode,Vector{Partition}}(),
)
    return depth_first_reconstruction(
        tree,
        conditioned_sample_node!,
        model,
        partition_list = partition_list,
        node_message_dict = node_message_dict,
    )
end
