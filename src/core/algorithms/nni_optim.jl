function nni_optim!(
        temp_message::Vector{Partition},
        message_to_set::Vector{Partition},
        node::GeneralFelNode,
        models,
        partition_list;
        acc_rule = (x,y) -> x > y)
    
    model_list = models(node)
    
    if isleafnode(node)
        return
    end

    #This bit of code should be identical to the regular downward pass...
    #-------------------
    
    for part in partition_list
        forward!(temp_message[part],node.parent_message[part],model_list[part],node)
    end
    @assert length(node.children) <= 2
    for i in 1:length(node.children)
        new_temp = deepcopy(temp_message) #Need to think of how to avoid this allocation. Same as in felsenstein_down
        sib_inds = sibling_inds(node.children[i])
        for part in partition_list
            combine!((node.children[i]).parent_message[part],[mess[part] for mess in node.child_messages[sib_inds]],true)
            combine!((node.children[i]).parent_message[part],[temp_message[part]],false)
        end
        #But calling branchlength_optim recursively...
        nni_optim!(new_temp,node.child_messages[i],node.children[i],models,partition_list; acc_rule = acc_rule)
    end
    #Then combine node.child_messages into node.message...
    for part in partition_list
        combine!(node.message[part],[mess[part] for mess in node.child_messages],true)
    end
    
    #But now we need to optimize the current node, and then prop back up to set your parents children message correctly.
    #-------------------
    if !isroot(node)
        nnid, exceed_sib, exceed_child = do_nni(node, temp_message, models; partition_list = partition_list, acc_rule = acc_rule)
        for part in partition_list
            combine!(node.message[part],[mess[part] for mess in node.child_messages],true)
            backward!(message_to_set[part],node.message[part],model_list[part],node)
            combine!(node.parent.message[part],[mess[part] for mess in node.parent.child_messages],true)
        end
    end
end

#Unsure if this is the best choice to handle the model,models, and model_func stuff.
function nni_optim!(
    temp_message::Vector{Partition},
    message_to_set::Vector{Partition},
    node::GeneralFelNode,
    models::Vector{<:BranchModel},
    partition_list;
    acc_rule = (x,y) -> x > y)
    nni_optim!(temp_message,message_to_set,node,x -> models,partition_list, acc_rule = acc_rule)
end
function nni_optim!(
    temp_message::Vector{Partition},
    message_to_set::Vector{Partition},
    node::GeneralFelNode,
    model::BranchModel,
    partition_list;
    acc_rule = (x,y) -> x > y)
    nni_optim!(temp_message,message_to_set,node,x -> [model],partition_list, acc_rule = acc_rule)
end

function do_nni(node, temp_message, models::F; 
        partition_list=1:length(node.message), acc_rule = (x,y) -> x > y) where F <: Function
    if length(node.children) == 0 || node.parent === nothing
        return false
    else
        temp_message2 = deepcopy(temp_message)
        model_list = models(node)
        #current score
        for part in partition_list
            backward!(temp_message[part],node.message[part],model_list[part],node)
            combine!(temp_message[part],[node.parent_message[part]],false)
        end
        #@toggleable_function assert_message_consistency(node, models, p = 0.01)

        curr_LL = sum([
                    total_LL(temp_message[part]) #+ 
                    #total_LL(node.message[part]) +
                    #total_LL(node.parent_message[part]) 
                for part in partition_list])

        max_LL = -Inf
        exceeded, exceed_sib, exceed_child = (false, 0, 0)

        for sib_ind in [x for x in 1:length(node.parent.children) if node.parent.children[x] != node]
            switch_LL = 0.0
            for child_ind in 1:length(node.children)
                for part in partition_list
                    #move the sibling message, after upward propogation, to temp_message to work with it
                    combine!(temp_message[part], [node.parent.child_messages[sib_ind][part]], true)

                    #combine this message, with all child messages of node except the index replaced
                    combine!(temp_message[part], [mess[part] for (i,mess) in enumerate(node.child_messages) if i != child_ind], false)

                    #prop up the message on the node up to its parent
                    backward!(temp_message2[part],temp_message[part],model_list[part],node)

                    #combine the message of the moved child
                    combine!(temp_message2[part], [node.child_messages[child_ind][part]], false)

                    #we now have both parts of message, propogated to the parent of node
                    #propogate it up one more step, then merge it with parent_message of parent
                    backward!(temp_message[part], temp_message2[part], model_list[part], node.parent)
                    combine!(temp_message[part], [node.parent.parent_message[part]], false)
                end

                switch_LL = sum([
                            total_LL(temp_message[part])
                        for part in partition_list])

                                            
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

            node.parent.child_messages[exceed_sib], node.child_messages[exceed_child] = node.child_messages[exceed_child], node.parent.child_messages[exceed_sib]

            return true, exceed_sib, exceed_child
        end
    end
end

function nni_optim!(tree::GeneralFelNode, models; partition_list = nothing,acc_rule = (x,y) -> x > y)
        temp_message = deepcopy(tree.message)
        message_to_set = deepcopy(tree.message)

        if partition_list === nothing
            partition_list = 1:length(tree.message)
        end

        nni_optim!(temp_message,message_to_set,tree,models,partition_list,acc_rule = acc_rule)
end
