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

#I need to add a version of this that takes a generic optimizer function and uses that instead of golden_section_maximize on just the branchlength.
#This is for cases where the user stores node-level parameters and wants to optimize them.
function branchlength_optim!(
    temp_messages::Vector{Vector{T}},
    tree::FelNode,
    models,
    partition_list,
    tol;
    bl_optimizer::UnivariateOpt = GoldenSectionOpt(),
    traversal = Iterators.reverse
) where {T <: Partition}

    stack = [(pop!(temp_messages), tree, 1, 1, true, true)]
    while !isempty(stack)
        temp_message, node, ind, lastind, first, down = pop!(stack)
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
                    #Temp must be constant between iterations for a node during down...
                    child_iter = traversal(1:length(node.children))
                    lastind = Base.first(child_iter) #(which is why we track the last child to be visited during down)
                    push!(stack, (Vector{T}(), node, ind, lastind, false, false)) #... but not up
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
                    #But calling branchlength_optim! recursively... (the iterative equivalent)
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
                    push!(temp_messages, temp_message)
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
            push!(temp_messages, temp_message)
        end
    end
end

#BM: Check if running felsenstein_down! makes a difference.
"""
    branchlength_optim!(tree::FelNode, models;  <keyword arguments>)

Uses golden section search, or optionally Brent's method, to optimize all branches recursively, maintaining the integrity of the messages.
Requires felsenstein!() to have been run first.
models can either be a single model (if the messages on the tree contain just one Partition) or an array of models, if the messages have >1 Partition, or 
a function that takes a node, and returns a Vector{<:BranchModel} if you need the models to vary from one branch to another.

# Keyword Arguments
- `partition_list=nothing`: (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over (but you probably want to optimize branch lengths with all models, the default option).
- `tol=1e-5`: absolute tolerance for the `bl_optimizer`.
- `bl_optimizer=GoldenSectionOpt()`: univariate branchlength optimizer, has Brent's method as an option by setting bl_optimizer=BrentsMethodOpt().
- `sort_tree=false`: determines if a [`lazysort!`](@ref) will be performed, which can reduce the amount of temporary messages that has to be initialized.
- `traversal=Iterators.reverse`: a function that determines the traversal, permutes an iterable.
- `shuffle=false`: do a randomly shuffled traversal, overrides `traversal`.
"""
function branchlength_optim!(tree::FelNode, models; partition_list = nothing, tol = 1e-5, bl_optimizer::UnivariateOpt = GoldenSectionOpt(), sort_tree = false, traversal = Iterators.reverse, shuffle = false)
    sort_tree && lazysort!(tree) #A lazysorted tree minimizes the amount of temp_messages needed
    temp_messages = [copy_message(tree.message)]

    if partition_list === nothing
        partition_list = 1:length(tree.message)
    end

    branchlength_optim!(temp_messages, tree, models, partition_list, tol, bl_optimizer=bl_optimizer, traversal = shuffle ? x -> sample(x, length(x), replace=false) : traversal)
end

#Overloading to allow for direct model and model vec inputs
branchlength_optim!(
    tree::FelNode,
    models::Vector{<:BranchModel};
    partition_list = nothing,
    tol = 1e-5,
    bl_optimizer::UnivariateOpt = GoldenSectionOpt(),
    sort_tree = false,
    traversal = Iterators.reverse,
    shuffle = false
) = branchlength_optim!(tree, x -> models, partition_list = partition_list, tol = tol, bl_optimizer = bl_optimizer, sort_tree = sort_tree, traversal = traversal, shuffle = shuffle)
branchlength_optim!(
    tree::FelNode,
    model::BranchModel;
    partition_list = nothing,
    tol = 1e-5,
    bl_optimizer::UnivariateOpt = GoldenSectionOpt(),
    sort_tree = false,
    traversal = Iterators.reverse,
    shuffle = false
) = branchlength_optim!(tree, x -> [model], partition_list = partition_list, tol = tol, bl_optimizer = bl_optimizer, sort_tree = sort_tree, traversal = traversal, shuffle = shuffle)
