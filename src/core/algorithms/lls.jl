
"""
total_LL(p::Partition)

If called on the root, it returns the log likelihood associated with that partition.
Can be overloaded for complex partitions without straightforward site log likelihoods.
"""
function total_LL(p::Partition)
    return sum(site_LLs(p))
end


function log_likelihood(tree::FelNode, model_func; partition_list = 1:length(tree.message))
    #Just combines root message with equilibrium message
    temp_message = copy_message(tree.parent_message[partition_list])
    models = model_func(tree)
    for (p, part) in enumerate(partition_list)
        #if tree.branchlength > 0.0 #Consider this - requires a (reasonable) assumption about model behavior.
        #Note: This reasonable assumption was violated by the SWM partition, which uses "forward!" to set the partition weights from the model.
            forward!(temp_message[p], tree.parent_message[part], models[part], tree)
        #end
        combine!(temp_message[p], tree.message[part])
    end
    return sum(total_LL.(temp_message))
end

"""
    log_likelihood(tree::FelNode, models; partition_list = nothing)

Computed the log likelihood of this tree. Requires felsenstein!() to have been run.
models can either be a single model (if the messages on the tree contain just one Partition) or an array of models, if the messages have >1 Partition, or 
a function that takes a node, and returns a Vector{<:BranchModel} if you need the models to vary from one branch to another.
partition_list (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over.
"""
log_likelihood(tree::FelNode, model::BranchModel; partition_list = 1:length(tree.message)) =
    log_likelihood(tree, x -> [model], partition_list = partition_list)
    
log_likelihood(
    tree::FelNode,
    models::Vector{<:BranchModel};
    partition_list = 1:length(tree.message),
) = log_likelihood(tree, x -> models, partition_list = partition_list)


"""
    log_likelihood!(tree::FelNode, models; partition_list = nothing)

First re-computes the upward felsenstein pass, and then computes the log likelihood of this tree.
models can either be a single model (if the messages on the tree contain just one Partition) or an array of models, if the messages have >1 Partition, or 
a function that takes a node, and returns a Vector{<:BranchModel} if you need the models to vary from one branch to another.
partition_list (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over.
"""
function log_likelihood!(tree::FelNode, models; partition_list = 1:length(tree.message))
    felsenstein!(tree, models, partition_list = partition_list)
    return log_likelihood(tree, models; partition_list = partition_list)
end
