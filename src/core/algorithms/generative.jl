
########Code related to simulation.


#Replaces an uncertain message with a sample from that distribution.
function sample_from_message!(message::Vector{<:Partition})
    sample_partition!.(message)
end

"""
sample_down!(node::GeneralFelNode,models,partition_list)

Generates samples under the model. The root.parent_message is taken as the starting distribution, and node.message contains the sampled messages.
"""
function sample_down!(node::GeneralFelNode,
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

"""
sample_down!(node::GeneralFelNode,models)

Generates samples under the model, matching the partition structure.
"""
function sample_down!(node::GeneralFelNode,
    models; partition_list = 1:length(node.message))
    sample_down!(node,models,partition_list)
end


"""
sample_down!(node::GeneralFelNode,models::Vector{<:BranchModel},partition_list)

Generates samples under the model, matching the partition structure.
"""
function sample_down!(node::GeneralFelNode,
    models::Vector{<:BranchModel},
    partition_list = 1:length(node.message))
sample_down!(node,x -> models,partition_list)
end

"""
sample_down!(node::GeneralFelNode,model::BranchModel,partition_list)

Generates samples under the model, matching the partition structure.
"""
function sample_down!(node::GeneralFelNode,
    models::BranchModel;
    partition_list = 1:length(node.message))
sample_down!(node,x -> [models],partition_list)
end