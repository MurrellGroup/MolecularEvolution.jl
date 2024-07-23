

function nni_optim_old!(
    temp_message::Vector{<:Partition},
    message_to_set::Vector{<:Partition},
    node::FelNode,
    models,
    partition_list;
    acc_rule = (x, y) -> x > y,
)    

    model_list = models(node)

    if isleafnode(node)
        return
    end

    #This bit of code should be identical to the regular downward pass...
    #-------------------

    for part in partition_list
        forward!(temp_message[part], node.parent_message[part], model_list[part], node)
    end
    @assert length(node.children) <= 2
    for i = 1:length(node.children)
        new_temp = copy_message(temp_message) #Need to think of how to avoid this allocation. Same as in felsenstein_down
        sib_inds = sibling_inds(node.children[i])
        for part in partition_list
            combine!(
                (node.children[i]).parent_message[part],
                [mess[part] for mess in node.child_messages[sib_inds]],
                true,
            )
            combine!((node.children[i]).parent_message[part], [temp_message[part]], false)
        end
        #But calling branchlength_optim recursively...
        nni_optim_old!(
            new_temp,
            node.child_messages[i],
            node.children[i],
            models,
            partition_list;
            acc_rule = acc_rule,
        )
    end
    #Then combine node.child_messages into node.message...
    for part in partition_list
        combine!(node.message[part], [mess[part] for mess in node.child_messages], true)
    end

    #But now we need to optimize the current node, and then prop back up to set your parents children message correctly.
    #-------------------
    if !isroot(node)
        nnid, exceed_sib, exceed_child = do_nni_old(
            node,
            temp_message,
            models;
            partition_list = partition_list,
            acc_rule = acc_rule,
        )
        for part in partition_list
            combine!(node.message[part], [mess[part] for mess in node.child_messages], true)
            backward!(message_to_set[part], node.message[part], model_list[part], node)
            combine!(
                node.parent.message[part],
                [mess[part] for mess in node.parent.child_messages],
                true,
            )
        end
    end
end

#Unsure if this is the best choice to handle the model,models, and model_func stuff.
function nni_optim_old!(
    temp_message::Vector{<:Partition},
    message_to_set::Vector{<:Partition},
    node::FelNode,
    models::Vector{<:BranchModel},
    partition_list;
    acc_rule = (x, y) -> x > y,
)
    nni_optim_old!(
        temp_message,
        message_to_set,
        node,
        x -> models,
        partition_list,
        acc_rule = acc_rule,
    )
end
function nni_optim_old!(
    temp_message::Vector{<:Partition},
    message_to_set::Vector{<:Partition},
    node::FelNode,
    model::BranchModel,
    partition_list;
    acc_rule = (x, y) -> x > y,
)
    nni_optim_old!(
        temp_message,
        message_to_set,
        node,
        x -> [model],
        partition_list,
        acc_rule = acc_rule,
    )
end

function do_nni_old(
    node,
    temp_message,
    models::F;
    partition_list = 1:length(node.message),
    acc_rule = (x, y) -> x > y,
) where {F<:Function}
    if length(node.children) == 0 || node.parent === nothing
        return false
    else
        temp_message2 = copy_message(temp_message)
        model_list = models(node)
        #current score
        for part in partition_list
            backward!(temp_message[part], node.message[part], model_list[part], node)
            combine!(temp_message[part], [node.parent_message[part]], false)
        end
        #@toggleable_function assert_message_consistency(node, models, p = 0.01)

        curr_LL = sum([total_LL(temp_message[part]) #+ 
        #total_LL(node.message[part]) +
        #total_LL(node.parent_message[part]) 
                       for part in partition_list])

        max_LL = -Inf
        exceeded, exceed_sib, exceed_child = (false, 0, 0)

        for sib_ind in
            [x for x in 1:length(node.parent.children) if node.parent.children[x] != node]
            switch_LL = 0.0
            for child_ind = 1:length(node.children)
                for part in partition_list
                    #move the sibling message, after upward propogation, to temp_message to work with it
                    combine!(
                        temp_message[part],
                        [node.parent.child_messages[sib_ind][part]],
                        true,
                    )

                    #combine this message, with all child messages of node except the index replaced
                    combine!(
                        temp_message[part],
                        [
                            mess[part] for
                            (i, mess) in enumerate(node.child_messages) if i != child_ind
                        ],
                        false,
                    )

                    #prop up the message on the node up to its parent
                    backward!(
                        temp_message2[part],
                        temp_message[part],
                        model_list[part],
                        node,
                    )

                    #combine the message of the moved child
                    combine!(
                        temp_message2[part],
                        [node.child_messages[child_ind][part]],
                        false,
                    )

                    #we now have both parts of message, propogated to the parent of node
                    #propogate it up one more step, then merge it with parent_message of parent
                    backward!(
                        temp_message[part],
                        temp_message2[part],
                        model_list[part],
                        node.parent,
                    )
                    combine!(temp_message[part], [node.parent.parent_message[part]], false)
                end

                switch_LL = sum([total_LL(temp_message[part]) for part in partition_list])


                if switch_LL > max_LL
                    exceed_sib = sib_ind
                    exceed_child = child_ind
                    max_LL = switch_LL
                end
            end
        end

        exceeded = acc_rule(max_LL, curr_LL)

        #do the actual move here, switching exceed child and exceed sib
        if !(exceeded)
            return false, exceed_sib, exceed_child
        else
            sib = node.parent.children[exceed_sib]
            child = node.children[exceed_child]

            child.parent = node.parent
            sib.parent = node

            node.children[exceed_child] = sib
            node.parent.children[exceed_sib] = child

            node.parent.child_messages[exceed_sib], node.child_messages[exceed_child] =
                node.child_messages[exceed_child], node.parent.child_messages[exceed_sib]

            return true, exceed_sib, exceed_child
        end
    end
end

"""
    nni_optim_old!(tree::FelNode, models; partition_list = nothing, tol = 1e-5)

Considers local branch swaps for all branches recursively, maintaining the integrity of the messages.
Requires felsenstein!() to have been run first.
models can either be a single model (if the messages on the tree contain just one Partition) or an array of models, if the messages have >1 Partition, or 
a function that takes a node, and returns a Vector{<:BranchModel} if you need the models to vary from one branch to another.
partition_list (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over (but you probably want to optimize tree topology with all models).
acc_rule allows you to specify a function that takes the current and proposed log likelihoods, and if true is returned the move is accepted.
"""
function nni_optim_old!(
    tree::FelNode,
    models;
    partition_list = nothing,
    acc_rule = (x, y) -> x > y,
)
    temp_message = copy_message(tree.message)
    message_to_set = copy_message(tree.message)

    if partition_list === nothing
        partition_list = 1:length(tree.message)
    end

    nni_optim_old!(
        temp_message,
        message_to_set,
        tree,
        models,
        partition_list,
        acc_rule = acc_rule,
    )
end

function nni_optim!(
    temp_message::Vector{<:Partition},
    message_to_set::Vector{<:Partition},
    node::FelNode,
    models,
    partition_list;
    acc_rule = (x, y) -> x > y,
    sampler = (x) -> (true, argmax(x[2:end]) + 1),
)

    model_list = models(node)

    if isleafnode(node)
        return
    end

    #This bit of code should be identical to the regular downward pass...
    #-------------------

    for part in partition_list
        forward!(temp_message[part], node.parent_message[part], model_list[part], node)
    end
    @assert length(node.children) <= 2
    for i = 1:length(node.children)
        new_temp = copy_message(temp_message) #Need to think of how to avoid this allocation. Same as in felsenstein_down
        sib_inds = sibling_inds(node.children[i])
        for part in partition_list
            combine!(
                (node.children[i]).parent_message[part],
                [mess[part] for mess in node.child_messages[sib_inds]],
                true,
            )
            combine!((node.children[i]).parent_message[part], [temp_message[part]], false)
        end
        #But calling branchlength_optim recursively...
        nni_optim!(
            new_temp,
            node.child_messages[i],
            node.children[i],
            models,
            partition_list;
            acc_rule = acc_rule,
            sampler = sampler,
        )
    end
    #Then combine node.child_messages into node.message...
    for part in partition_list
        combine!(node.message[part], [mess[part] for mess in node.child_messages], true)
    end

    #But now we need to optimize the current node, and then prop back up to set your parents children message correctly.
    #-------------------
    if !isroot(node)
        nnid, exceed_sib, exceed_child = do_nni(
            node,
            temp_message,
            models;
            partition_list = partition_list,
            acc_rule = acc_rule,
            sampler = sampler,
        )
        for part in partition_list
            combine!(node.message[part], [mess[part] for mess in node.child_messages], true)
            backward!(message_to_set[part], node.message[part], model_list[part], node)
            combine!(
                node.parent.message[part],
                [mess[part] for mess in node.parent.child_messages],
                true,
            )
        end
    end
end
#Unsure if this is the best choice to handle the model,models, and model_func stuff.
function nni_optim!(
    temp_message::Vector{<:Partition},
    message_to_set::Vector{<:Partition},
    node::FelNode,
    models::Vector{<:BranchModel},
    partition_list;
    acc_rule = (x, y) -> x > y,
    sampler = (x) ->  (true, argmax(x[2:end]) + 1),
)
    nni_optim!(
        temp_message,
        message_to_set,
        node,
        x -> models,
        partition_list,
        acc_rule = acc_rule,
        sampler = sampler,
    )
end
function nni_optim!(
    temp_message::Vector{<:Partition},
    message_to_set::Vector{<:Partition},
    node::FelNode,
    model::BranchModel,
    partition_list;
    acc_rule = (x, y) -> x > y,
    sampler = (x) ->  (true, argmax(x[2:end]) + 1),
)
    nni_optim!(
        temp_message,
        message_to_set,
        node,
        x -> [model],
        partition_list,
        acc_rule = acc_rule,
        sampler = sampler,
    )
end

function do_nni(
    node,
    temp_message,
    models::F;
    partition_list = 1:length(node.message),
    acc_rule = (x, y) -> x > y,
    sampler = (x) ->  (true, argmax(x[2:end]) + 1),
    ) where {F<:Function}

    if length(node.children) == 0 || node.parent === nothing
        return false
    else
        temp_message2 = copy_message(temp_message)
        model_list = models(node)
        #current score
        for part in partition_list
            backward!(temp_message[part], node.message[part], model_list[part], node)
            combine!(temp_message[part], [node.parent_message[part]], false)
        end
        #@toggleable_function assert_message_consistency(node, models, p = 0.01)

        curr_LL = sum([total_LL(temp_message[part]) #+ 
        #total_LL(node.message[part]) +
        #total_LL(node.parent_message[part]) 
                       for part in partition_list])
 
        changed = false
        nni_LLs = [curr_LL]
        nni_configs = [(0,0)]

        for sib_ind in
            [x for x in 1:length(node.parent.children) if node.parent.children[x] != node]
            switch_LL = 0.0
            for child_ind = 1:length(node.children)
                for part in partition_list
                    #move the sibling message, after upward propogation, to temp_message to work with it
                    combine!(
                        temp_message[part],
                        [node.parent.child_messages[sib_ind][part]],
                        true,
                    )

                    #combine this message, with all child messages of node except the index replaced
                    combine!(
                        temp_message[part],
                        [
                            mess[part] for
                            (i, mess) in enumerate(node.child_messages) if i != child_ind
                        ],
                        false,
                    )

                    #prop up the message on the node up to its parent
                    backward!(
                        temp_message2[part],
                        temp_message[part],
                        model_list[part],
                        node,
                    )

                    #combine the message of the moved child
                    combine!(
                        temp_message2[part],
                        [node.child_messages[child_ind][part]],
                        false,
                    )

                    #we now have both parts of message, propogated to the parent of node
                    #propogate it up one more step, then merge it with parent_message of parent
                    backward!(
                        temp_message[part],
                        temp_message2[part],
                        model_list[part],
                        node.parent,
                    )
                    combine!(temp_message[part], [node.parent.parent_message[part]], false)
                end

                
                LL = sum([total_LL(temp_message[part]) for part in partition_list])

                push!(nni_LLs, LL)
                push!(nni_configs, (sib_ind, child_ind))

            end
        end

        changed, sampled_config_ind = sampler(nni_LLs)
        sampled_config_LL = nni_LLs[sampled_config_ind]
        (sampled_sib_ind, sampled_child_ind) = nni_configs[sampled_config_ind]
        changed = acc_rule(sampled_config_LL, curr_LL) && changed
        

        # println("changed: ", changed, " inds, ", (sampled_sib_ind, sampled_child_ind), " sampled config ind, ", sampled_config_ind)
        # println(" nni_config, ", nni_configs)
        # println("softmax ", softmax(nni_LLs))

        #do the actual move here
        if !(changed)
            return false, sampled_sib_ind, sampled_child_ind
        else
            sib = node.parent.children[sampled_sib_ind]
            child = node.children[sampled_child_ind]

            child.parent = node.parent
            sib.parent = node

            node.children[sampled_child_ind] = sib
            node.parent.children[sampled_sib_ind] = child

            node.parent.child_messages[sampled_sib_ind], node.child_messages[sampled_child_ind] =
                node.child_messages[sampled_child_ind], node.parent.child_messages[sampled_sib_ind]

            return true, sampled_sib_ind, sampled_child_ind
        end
    end
end

"""
    nni_optim!(tree::FelNode, models; partition_list = nothing, tol = 1e-5)

Considers local branch swaps for all branches recursively, maintaining the integrity of the messages.
Requires felsenstein!() to have been run first.
models can either be a single model (if the messages on the tree contain just one Partition) or an array of models, if the messages have >1 Partition, or 
a function that takes a node, and returns a Vector{<:BranchModel} if you need the models to vary from one branch to another.
partition_list (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over (but you probably want to optimize tree topology with all models).
acc_rule allows you to specify a function that takes the current and proposed log likelihoods, and if true is returned the move is accepted.
sampler allows you to randomly select a NNI based on the vector of log likelihoods of the possible interchanges, note that the current log likelihood is at index 1. 
"""
function nni_optim!(
    tree::FelNode,
    models;
    partition_list = nothing,
    acc_rule = (x, y) -> x > y,
    sampler = (x) ->  (true, argmax(x[2:end]) + 1),
)
    temp_message = copy_message(tree.message)
    message_to_set = copy_message(tree.message)

    if partition_list === nothing
        partition_list = 1:length(tree.message)
    end

    nni_optim!(
        temp_message,
        message_to_set,
        tree,
        models,
        partition_list,
        acc_rule = acc_rule,
        sampler = sampler,
    )
end

#Overloading to allow for direct model and model vec inputs
nni_optim!(
    tree::FelNode,
    models::Vector{<:BranchModel};
    partition_list = nothing,
    acc_rule = (x, y) -> x > y,
    sampler = (x) ->  (true, argmax(x[2:end]) + 1),
) = nni_optim!(tree, x -> models, partition_list=partition_list, acc_rule=acc_rule, sampler = sampler)
nni_optim!(
    tree::FelNode,
    model::BranchModel;
    partition_list = nothing,
    acc_rule = (x, y) -> x > y,
    sampler = (x) ->  (true, argmax(x[2:end]) + 1),
) = nni_optim!(tree, x -> [model], partition_list=partition_list, acc_rule=acc_rule, sampler = sampler)


