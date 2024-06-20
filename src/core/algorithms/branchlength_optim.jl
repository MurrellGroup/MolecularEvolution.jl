#Model list should be a list of P matrices.
function branch_LL_up(
    bl::Real,
    temp_message::Vector{<:Partition},
    node::FelNode,
    model_list::Vector{<:BranchModel},
    partition_list,
)
    #Questionable:
    backup = node.branchlength
    node.branchlength = bl
    for part in partition_list
        backward!(temp_message[part], node.message[part], model_list[part], node)
        combine!(temp_message[part], [node.parent_message[part]], false)
    end
    tot_LL = sum([total_LL(temp_message[part]) for part in partition_list])
    #Debatable, but we're going to revert for now:
    node.branchlength = backup
    return tot_LL
end

#TODO: make use of a lazysort! and a message_stack (memoryblocks)
#TODO: use the exact same principle in nni_optim!
function branchlength_optim_iter!(
    temp_message::Vector{<:Partition},
    tree::FelNode,
    models,
    partition_list,
    tol;
    bl_optimizer::UnivariateOpt = GoldenSectionOpt()
)

    stack = [(temp_message, tree, 1, true, true)]
    while !isempty(stack)
        temp_message, node, ind, first, down = pop!(stack)
        #We start out with a regular downward pass...
        #(except for some extra bookkeeping to track if node is visited for the first time)
        #-------------------
        if !isleafnode(node)
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
                    push!(stack, (temp_message, node, ind, false, false))
                    for i = Iterators.reverse(1:length(node.children)) #Iterative reverse <=> Recursive non-reverse, also optimal for lazysort!??
                        push!(stack, (temp_message, node, i, false, true))
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
                    #But calling branchlength_optim! recursively... (the iterative equivalent)
                    push!(stack, (copy_message(temp_message), node.children[ind], ind, true, true))
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
                    model_list = models(node)
                    fun = x -> branch_LL_up(x, temp_message, node, model_list, partition_list)
                    opt = univariate_maximize(fun, 0 + tol, 1 - tol, unit_transform, bl_optimizer, tol)
                    if fun(opt) > fun(node.branchlength)
                        node.branchlength = opt
                    end
                    #Consider checking for improvement, and bailing if none.
                    #Then we need to set the "message_to_set", which is node.parent.child_messages[but_the_right_one]
                    for part in partition_list
                        backward!(
                            node.parent.child_messages[ind][part],
                            node.message[part],
                            model_list[part],
                            node,
                        )
                    end
                end
            end
        else
            #But now we need to optimize the current node, and then prop back up to set your parents children message correctly.
            #-------------------
            model_list = models(node)
            fun = x -> branch_LL_up(x, temp_message, node, model_list, partition_list)
            opt = univariate_maximize(fun, 0 + tol, 1 - tol, unit_transform, bl_optimizer, tol)
            if fun(opt) > fun(node.branchlength)
                node.branchlength = opt
            end
            #Consider checking for improvement, and bailing if none.
            #Then we need to set the "message_to_set", which is node.parent.child_messages[but_the_right_one]
            for part in partition_list
                backward!(
                    node.parent.child_messages[ind][part],
                    node.message[part],
                    model_list[part],
                    node,
                )
            end
        end
    end
end

#I need to add a version of this that takes a generic optimizer function and uses that instead of golden_section_maximize on just the branchlength.
#This is for cases where the user stores node-level parameters and wants to optimize them.
function branchlength_optim!(
    temp_message::Vector{<:Partition},
    message_to_set::Vector{<:Partition},
    node::FelNode,
    models,
    partition_list,
    tol;
    bl_optimizer::UnivariateOpt = GoldenSectionOpt()
)

    #This bit of code should be identical to the regular downward pass...
    #-------------------
    if !isleafnode(node)
        model_list = models(node)
        for part in partition_list
            forward!(temp_message[part], node.parent_message[part], model_list[part], node)
        end
        for i = 1:length(node.children)
            new_temp = copy_message(temp_message) #Need to think of how to avoid this allocation. Same as in felsenstein_down
            sib_inds = sibling_inds(node.children[i])
            for part in partition_list
                combine!(
                    (node.children[i]).parent_message[part],
                    [mess[part] for mess in node.child_messages[sib_inds]],
                    true,
                )
                combine!(
                    (node.children[i]).parent_message[part],
                    [temp_message[part]],
                    false,
                )
            end
            #But calling branchlength_optim recursively...
            branchlength_optim!(
                new_temp,
                node.child_messages[i],
                node.children[i],
                models,
                partition_list,
                tol,
                bl_optimizer=bl_optimizer
            )
        end
        #Then combine node.child_messages into node.message...
        for part in partition_list
            combine!(node.message[part], [mess[part] for mess in node.child_messages], true)
        end
    end
    #But now we need to optimize the current node, and then prop back up to set your parents children message correctly.
    #-------------------
    if !isroot(node)
        model_list = models(node)
        fun = x -> branch_LL_up(x, temp_message, node, model_list, partition_list)
        opt = univariate_maximize(fun, 0 + tol, 1 - tol, unit_transform, bl_optimizer, tol)
        if fun(opt) > fun(node.branchlength)
            node.branchlength = opt
        end
        #Consider checking for improvement, and bailing if none.
        #Then we need to set the "message_to_set", which is node.parent.child_messages[but_the_right_one]
        for part in partition_list
            backward!(message_to_set[part], node.message[part], model_list[part], node)
        end
    end
    #For debugging:
    #println("$(node.nodeindex):$(node.branchlength)")
end

function branchlength_optim_iter!(tree::FelNode, models; partition_list = nothing, tol = 1e-5, bl_optimizer::UnivariateOpt = GoldenSectionOpt())
    temp_message = copy_message(tree.message)

    if partition_list === nothing
        partition_list = 1:length(tree.message)
    end

    branchlength_optim_iter!(temp_message, tree, models, partition_list, tol, bl_optimizer=bl_optimizer)
end

#BM: Check if running felsenstein_down! makes a difference.
"""
    branchlength_optim!(tree::FelNode, models; partition_list = nothing, tol = 1e-5, bl_optimizer::UnivariateOpt = GoldenSectionOpt())

Uses golden section search, or optionally Brent's method, to optimize all branches recursively, maintaining the integrity of the messages.
Requires felsenstein!() to have been run first.
models can either be a single model (if the messages on the tree contain just one Partition) or an array of models, if the messages have >1 Partition, or 
a function that takes a node, and returns a Vector{<:BranchModel} if you need the models to vary from one branch to another.
partition_list (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over (but you probably want to optimize branch lengths with all models).
tol is the absolute tolerance for the bl_optimizer which defaults to golden section search, and has Brent's method as an option by setting bl_optimizer=BrentsMethodOpt().
"""
function branchlength_optim!(tree::FelNode, models; partition_list = nothing, tol = 1e-5, bl_optimizer::UnivariateOpt = GoldenSectionOpt())
    temp_message = copy_message(tree.message)
    message_to_set = copy_message(tree.message)

    if partition_list === nothing
        partition_list = 1:length(tree.message)
    end

    branchlength_optim!(temp_message, message_to_set, tree, models, partition_list, tol, bl_optimizer=bl_optimizer)
end

#Overloading to allow for direct model and model vec inputs
branchlength_optim!(
    tree::FelNode,
    models::Vector{<:BranchModel};
    partition_list = nothing,
    tol = 1e-5,
    bl_optimizer::UnivariateOpt = GoldenSectionOpt()
) = branchlength_optim!(tree, x -> models, partition_list = partition_list, tol = tol, bl_optimizer=bl_optimizer)
branchlength_optim!(
    tree::FelNode,
    model::BranchModel;
    partition_list = nothing,
    tol = 1e-5,
    bl_optimizer::UnivariateOpt = GoldenSectionOpt()
) = branchlength_optim!(tree, x -> [model], partition_list = partition_list, tol = tol, bl_optimizer=bl_optimizer)
