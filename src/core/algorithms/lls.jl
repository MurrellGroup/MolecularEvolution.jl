
function total_LL(p::Partition)
    return sum(site_LLs(p))
end

"""
    log_likelihood(tree::FelNode; partition_list = nothing)

Computed the log likelihood of this tree. Requires felsenstein!() to have been run.
"""
function log_likelihood(tree::FelNode, model_func; partition_list = 1:length(tree.message))
    #Just combines root message with equilibrium message
    temp_message = deepcopy(tree.parent_message[partition_list])
    models = model_func(tree)
    for (p,part) in enumerate(partition_list)
        if tree.branchlength > 0.0 #Consider this - requires a (reasonable) assumption about model behavior.
            forward!(temp_message[p], tree.parent_message[part], models[part],tree)
        end
        #BM: I think lots of code will break if eg. partition list = [2]
        combine!(temp_message[p], tree.message[part])
    end
    return sum(total_LL.(temp_message))
end

log_likelihood(tree::FelNode, model::BranchModel; partition_list = 1:length(tree.message)) = log_likelihood(tree, x -> [model], partition_list = partition_list)
log_likelihood(tree::FelNode, models::Vector{<:BranchModel}; partition_list = 1:length(tree.message)) = log_likelihood(tree, x -> models, partition_list = partition_list)


"""
    log_likelihood!(tree::FelNode, models; partition_list = nothing)

First re-computes the upward felsenstein pass, and then computes the log likelihood of this tree.
"""
function log_likelihood!(tree::FelNode, models; partition_list = 1:length(tree.message))
    felsenstein!(tree, models, partition_list = partition_list)
    return log_likelihood(tree, models; partition_list=partition_list)
end
