#Model list should be a list of P matrices.
function branch_LL_up(bl::Real,
        temp_message::Vector{Partition},
        node::FelNode,
        model_list::Vector{<:BranchModel},
        partition_list)
    #Questionable:
    backup = node.branchlength
    node.branchlength = bl
    for part in partition_list
        backward!(temp_message[part],node.message[part],model_list[part],node)
        combine!(temp_message[part],[node.parent_message[part]],false)
    end
    tot_LL = sum([total_LL(temp_message[part]) for part in partition_list])
    #Debatable, but we're going to revert for now:
    node.branchlength = backup
    return tot_LL
end

#I need to add a version of this that takes a generic optimizer function and uses that instead of golden_section_maximize on just the branchlength.
#This is for cases where the user stores node-level parameters and wants to optimize them.
function branchlength_optim!(
        temp_message::Vector{Partition},
        message_to_set::Vector{Partition},
        node::FelNode,
        models,
        partition_list,
        tol)

    #This bit of code should be identical to the regular downward pass...
    #-------------------
    if !isleafnode(node)
        model_list = models(node)
        for part in partition_list
            forward!(temp_message[part],node.parent_message[part],model_list[part],node)
        end
        for i in 1:length(node.children)
            new_temp = deepcopy(temp_message) #Need to think of how to avoid this allocation. Same as in felsenstein_down
            sib_inds = sibling_inds(node.children[i])
            for part in partition_list
                combine!((node.children[i]).parent_message[part],[mess[part] for mess in node.child_messages[sib_inds]],true)
                combine!((node.children[i]).parent_message[part],[temp_message[part]],false)
            end
            #But calling branchlength_optim recursively...
            branchlength_optim!(new_temp,node.child_messages[i],node.children[i],models,partition_list,tol)
        end
        #Then combine node.child_messages into node.message...
        for part in partition_list
            combine!(node.message[part],[mess[part] for mess in node.child_messages],true)
        end
    end
    #But now we need to optimize the current node, and then prop back up to set your parents children message correctly.
    #-------------------
    if !isroot(node)
        model_list = models(node)
        fun = x->branch_LL_up(x,temp_message,node,model_list,partition_list)
        opt = golden_section_maximize(fun, 0+tol, 1-tol,unit_transform,tol)
        if fun(opt) > fun(node.branchlength)
            node.branchlength = opt
        end
        #Consider checking for improvement, and bailing if none.
        #Then we need to set the "message_to_set", which is node.parent.child_messages[but_the_right_one]
        for part in partition_list
            backward!(message_to_set[part],node.message[part],model_list[part],node)
        end
    end
    #For debugging:
    #println("$(node.nodeindex):$(node.branchlength)")
end


"""
    branchlength_optim!(tree::FelNode, models; partition_list = nothing, tol = 1e-5)

Uses golden section search to optimize all branches recursively, maintaining the integrity of the messages.
Requires felsenstein!() to have been run first.
"""
function branchlength_optim!(tree::FelNode, models; partition_list = nothing, tol = 1e-5)
        temp_message = deepcopy(tree.message)
        message_to_set = deepcopy(tree.message)

        if partition_list === nothing
            partition_list = 1:length(tree.message)
        end

        branchlength_optim!(temp_message,message_to_set,tree,models,partition_list,tol)
end

#Overloading to allow for direct model and model vec inputs
branchlength_optim!(tree::FelNode, models::Vector{<:BranchModel} ; partition_list = nothing, tol = 1e-5) = branchlength_optim!(tree, x -> models, partition_list = partition_list, tol = tol)
branchlength_optim!(tree::FelNode, model::BranchModel; partition_list = nothing, tol = 1e-5) = branchlength_optim!(tree, x -> [model], partition_list = partition_list, tol = tol)