#After a do_nni, we have to update parent_message if we want to continue down (assume that temp_message is the forwarded parent.parent_message)
function update_parent_message!(
    node::FelNode,
    temp_message::Vector{<:Partition};
    partition_list = 1:length(node.message),
)
    sib_inds = sibling_inds(node)
    for part in partition_list
        combine!(
            node.parent_message[part],
            [mess[part] for mess in node.parent.child_messages[sib_inds]],
            true,
        )
        combine!(
            node.parent_message[part],
            [temp_message[part]],
            false,
        )
    end
end

function nni_optim!(
    temp_messages::Vector{Vector{T}},
    tree::FelNode,
    models,
    partition_list;
    selection_rule = x -> argmax(x),
    traversal = Iterators.reverse
) where {T <: Partition}

    #Consider a NamedTuple/struct
    stack = [(pop!(temp_messages), tree, 1, 1, true, true)]
    while !isempty(stack)
        temp_message, node, ind, lastind, first, down = pop!(stack)
        #We start out with a regular downward pass...
        #(except for some extra bookkeeping to track if node is visited for the first time)
        #-------------------
        if isleafnode(node)
            push!(temp_messages, temp_message)
            continue
        end
        if down
            if first
                model_list = models(node)
                for part in partition_list
                    forward!(
                        temp_message[part],
                        node.parent_message[part],
                        model_list[part],
                        node,
                    )
                end
                @assert length(node.children) <= 2
                #Temp must be constant between iterations for a node during down...
                child_iter = traversal(1:length(node.children))
                lastind = Base.first(child_iter) #(which is why we track the last child to be visited during down)
                push!(stack, (Vector{T}(), node, ind, lastind, true, false)) #... but not up
                for i = child_iter #Iterative reverse <=> Recursive non-reverse, also optimal for lazysort!??
                    push!(stack, (temp_message, node, i, lastind, false, true))
                end
            end
            if !first
                sib_inds = sibling_inds(node.children[ind])
                for part in partition_list
                    combine!(
                        (node.children[ind]).parent_message[part],
                        [mess[part] for mess in node.child_messages[sib_inds]],
                        true,
                    )
                    combine!(
                        (node.children[ind]).parent_message[part],
                        [temp_message[part]],
                        false,
                    )
                end
                #But calling nni_optim! recursively... (the iterative equivalent)
                push!(stack, (safepop!(temp_messages, temp_message), node.children[ind], ind, lastind, true, true)) #first + down combination => safepop!
                ind == lastind && push!(temp_messages, temp_message) #We no longer need constant temp
            end
        end
        if !down
            #Then combine node.child_messages into node.message...
            for part in partition_list
                combine!(node.message[part], [mess[part] for mess in node.child_messages], true)
            end
            #But now we need to optimize the current node, and then prop back up to set your parents children message correctly.
            #-------------------
            if !isroot(node)
                temp_message = pop!(temp_messages)
                model_list = models(node)
                if first #We only do_nni first up
                    nnid, sampled_sib_ind, sampled_child_ind = do_nni(
                        node,
                        temp_message,
                        models;
                        partition_list = partition_list,
                        selection_rule = selection_rule,
                    )
                    if nnid && last(last(stack)) #We nnid a sibling that hasn't been visited (then, down would be true in the next iter)...
                        #... and now we want to continue down the nnid sibling (now a child to node)
                        push!(temp_messages, temp_message)
                        temp_message = Base.first(last(stack)) #The forwarded parent.parent_message
                        #First we update the parent_message...
                        update_parent_message!(
                            node,
                            temp_message;
                            partition_list = partition_list,
                        )
                        #... then we forward the updated parent_message (this resembles a first down)
                        model_list = models(node)
                        for part in partition_list
                            forward!(
                                temp_message[part],
                                node.parent_message[part],
                                model_list[part],
                                node,
                            )
                        end
                        pop!(stack)
                        push!(stack, (Vector{T}(), node, ind, lastind, false, false)) #When we're going up a second time, we no longer need a temp
                        push!(stack, (temp_message, node, sampled_child_ind, sampled_child_ind, false, true)) #Go to the "new" child - the "new" lastind
                        continue #Don't fel-up yet
                    end
                end
                for part in partition_list
                    combine!(node.message[part], [mess[part] for mess in node.child_messages], true)
                    backward!(node.parent.child_messages[ind][part], node.message[part], model_list[part], node)
                    combine!(
                        node.parent.message[part],
                        [mess[part] for mess in node.parent.child_messages],
                        true,
                    )
                end
                push!(temp_messages, temp_message)
            end
        end
    end
end

#Unsure if this is the best choice to handle the model,models, and model_func stuff.
function nni_optim!(
    temp_messages::Vector{Vector{T}},
    tree::FelNode,
    models::Vector{<:BranchModel},
    partition_list;
    selection_rule = x -> argmax(x),
    traversal = Iterators.reverse,
) where {T <: Partition}
    nni_optim!(
        temp_messages,
        tree,
        x -> models,
        partition_list,
        selection_rule = selection_rule,
        traversal = traversal,
    )
end
function nni_optim!(
    temp_messages::Vector{Vector{T}},
    tree::FelNode,
    model::BranchModel,
    partition_list;
    selection_rule = x -> argmax(x),
    traversal = Iterators.reverse,

) where {T <: Partition}
    nni_optim!(
        temp_messages,
        tree,
        x -> [model],
        partition_list,
        selection_rule = selection_rule,
        traversal = traversal,
    )
end

function do_nni(
    node,
    temp_message,
    models::F;
    partition_list = 1:length(node.message),
    selection_rule = x -> argmax(x),
) where {F<:Function}
    if length(node.children) == 0 || node.parent === nothing
        return false
    else
        temp_message2 = copy_message(temp_message) #Make use of temp_messages here
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

        change = false
        nni_LLs = [curr_LL]
        nni_configs = [(0,0)]


        

        for sib_ind in
            [x for x in 1:length(node.parent.children) if node.parent.children[x] != node]

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

        sampled_config_ind = selection_rule(nni_LLs)
        change = sampled_config_ind != 1
        (sampled_sib_ind, sampled_child_ind) = nni_configs[sampled_config_ind]

        #do the actual move here, switching sampled_child_in and sampled_sib_ind
        if !(change)
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
    nni_optim!(tree::FelNode, models; <keyword arguments>)

Considers local branch swaps for all branches recursively, maintaining the integrity of the messages.
Requires felsenstein!() to have been run first.
models can either be a single model (if the messages on the tree contain just one Partition) or an array of models, if the messages have >1 Partition, or 
a function that takes a node, and returns a Vector{<:BranchModel} if you need the models to vary from one branch to another.

# Keyword Arguments
- `partition_list=nothing`: (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over (but you probably want to optimize tree topology with all models, the default option).
- `selection_rule = x -> argmax(x)`: a function that takes the current and proposed log likelihoods and selects a nni configuration. Note that the current log likelihood is stored at x[1].
- `sort_tree=false`: determines if a [`lazysort!`](@ref) will be performed, which can reduce the amount of temporary messages that has to be initialized.
- `traversal=Iterators.reverse`: a function that determines the traversal, permutes an iterable.
- `shuffle=false`: do a randomly shuffled traversal, overrides `traversal`.
"""
function nni_optim!(
    tree::FelNode,
    models;
    partition_list = nothing,
    selection_rule = x -> argmax(x),
    sort_tree = false,
    traversal = Iterators.reverse,
    shuffle = false
)
    sort_tree && lazysort!(tree) #A lazysorted tree minimizes the amount of temp_messages needed
    temp_messages = [copy_message(tree.message)]

    if partition_list === nothing
        partition_list = 1:length(tree.message)
    end

    nni_optim!(
        temp_messages,
        tree,
        models,
        partition_list,
        selection_rule = selection_rule,
        traversal = shuffle ? x -> sample(x, length(x), replace=false) : traversal
    )
end