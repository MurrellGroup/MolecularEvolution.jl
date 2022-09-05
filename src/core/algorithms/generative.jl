
########Code related to simulation.


"""
    sample_from_message!(message::Vector{<:Partition})

#Replaces an uncertain message with a sample from the distribution represented by each partition.
"""
function sample_from_message!(message::Vector{<:Partition})
    sample_partition!.(message)
end

"""
sample_down!(root::FelNode,models,partition_list)

Generates samples under the model. The root.parent_message is taken as the starting distribution, and node.message contains the sampled messages.
models can either be a single model (if the messages on the tree contain just one Partition) or an array of models, if the messages have >1 Partition, or 
a function that takes a node, and returns a Vector{<:BranchModel} if you need the models to vary from one branch to another.
partition_list (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over.
"""
function sample_down!(node::FelNode,
    models,
    partition_list)
    model_list = models(node)
    for part in partition_list
        if isroot(node)
            forward!(node.message[part],node.parent_message[part],model_list[part],node)
        else
            forward!(node.message[part],node.parent.message[part],model_list[part],node)
        end
        sample_partition!(node.message[part])
    end
    if !isleafnode(node)
        for child in node.children
            sample_down!(child,models,partition_list)
        end
    end
end

function sample_down!(node::FelNode,
    models; partition_list = 1:length(node.message))
    sample_down!(node,models,partition_list)
end

function sample_down!(node::FelNode,
    models::Vector{<:BranchModel},
    partition_list = 1:length(node.message))
    sample_down!(node,x -> models,partition_list)
end

function sample_down!(node::FelNode,
    models::BranchModel;
    partition_list = 1:length(node.message))
    sample_down!(node,x -> [models],partition_list)
end