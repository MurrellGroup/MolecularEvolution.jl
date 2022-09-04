"""
    felsenstein!(node::GeneralFelNode, models; partition_list = nothing)

Should usually be called on the root of the tree. Propagates Felsenstein pass up from the tips to the root.
models must be a function that takes a node, and returns a Vector{BranchModel}. This lets you control which models gets used for which branch.
partition_list (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over.
"""
function felsenstein!(tree::GeneralFelNode, models; partition_list = 1:length(tree.message))
    stack = [(tree,1, true)]
    #Note to future self: I tried replacing this with an actual stack, but it wasn't a big enough perf diff to justify the extra dependency
    #stack = Stack{Tuple{GeneralFelNode, Int64, Bool}}()
    #push!(stack, (tree,1, true))
    node = nothing
    while length(stack) > 0
        node, ind, first = pop!(stack)
        #println(node.nodeindex)
        if !isleafnode(node)
            if first
                push!(stack, (node,ind,false))
                for i in 1:length(node.children)
                    push!(stack, (node.children[i],i,true))
                end
            end
            if !first
                for part in partition_list
                    #Combine child messages into node message.
                    combine!(node.message[part],[mess[part] for mess in node.child_messages],true)
                end
                if !isroot(node)
                    #Get the model list for the current branch.
                    model_list = models(node)
                    for part in partition_list
                        backward!(node.parent.child_messages[ind][part],node.message[part],model_list[part],node)
                    end
                end
            end
        else
            model_list = models(node)
            for part in partition_list
                backward!(node.parent.child_messages[ind][part],node.message[part],model_list[part],node)
            end
        end
    end
end

"""
    felsenstein!(node::GeneralFelNode, models::Vector{<:BranchModel}; partition_list = nothing)

Should usually be called on the root of the tree. Propagates Felsenstein pass up from the tips to the root.
models can be a single model, or a vector of models, or a function that associates a node with a model.
partition_list (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over.
"""
function felsenstein!(tree::GeneralFelNode, models::Vector{<:BranchModel}; partition_list = 1:length(tree.message))
    felsenstein!(tree, x -> models, partition_list = partition_list)
end

"""
    felsenstein!(node::GeneralFelNode, model::BranchModel; partition_list = nothing)

Should usually be called on the root of the tree. Propagates Felsenstein pass up from the tips to the root.
models can be a single model, or a vector of models, or a function that associates a node with a model.
partition_list (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over.
"""
function felsenstein!(tree::GeneralFelNode, model::BranchModel; partition_list = 1:length(tree.message))
    felsenstein!(tree, x -> [model], partition_list = partition_list)
end

"""
    felsenstein_down!(node::GeneralFelNode, models; partition_list = 1:length(tree.message), temp_message = deepcopy(tree.message))

Should usually be called on the root of the tree. Propagates Felsenstein pass down from the root to the tips.
felsenstein!() should usually be called first.
models must be a function that takes a node, and returns a Vector{BranchModel}. This lets you control which models gets used for which branch.
partition_list (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over.
"""
function felsenstein_down!(tree::GeneralFelNode, models; partition_list = 1:length(tree.message), temp_message = deepcopy(tree.message))
    stack = [tree]
	#curr = nothing
	while length(stack) > 0
		node = pop!(stack)
		model_list = models(node)

		if !isleafnode(node)
        	for part in partition_list
        	    forward!(temp_message[part],node.parent_message[part],model_list[part],node)
        	end
			#Note that this uses a nested loop over children and siblings.
			#Thats fine for binary trees, but will blow up when we have massive polytomies.
			#Which is stupid, because those are equivalent to binary trees with zero branch lengths.
			#Either make this cleverer, or always binarize the trees.
			for i in 1:length(node.children)
				sib_inds = sibling_inds(node.children[i])
				for part in partition_list
					if sum(sib_inds) > 0
						combine!((node.children[i]).parent_message[part],[mess[part] for mess in node.child_messages[sib_inds]],true)
						combine!((node.children[i]).parent_message[part],[temp_message[part]],false)
					else
						combine!((node.children[i]).parent_message[part],[temp_message[part]],true)
					end
				end
				push!(stack, node.children[i])
			end
		end
	end
end

"""
    felsenstein_down!(node::GeneralFelNode, models; partition_list = 1:length(tree.message), temp_message = deepcopy(tree.message))

Should usually be called on the root of the tree. Propagates Felsenstein pass down from the root to the tips.
felsenstein!() should usually be called first.
models can be a single model, or a vector of models, or a function that associates a node with a model.
partition_list (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over.
"""
function felsenstein_down!(tree::GeneralFelNode, models::Vector{<:BranchModel}; partition_list = 1:length(tree.message), temp_message = deepcopy(tree.message))
    felsenstein_down!(tree, x -> models, partition_list = partition_list, temp_message = temp_message)
end

"""
    felsenstein_down!(node::GeneralFelNode, models; partition_list = 1:length(tree.message), temp_message = deepcopy(tree.message))

Should usually be called on the root of the tree. Propagates Felsenstein pass down from the root to the tips.
felsenstein!() should usually be called first.
models can be a single model, or a vector of models, or a function that associates a node with a model.
partition_list (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over.
"""
function felsenstein_down!(tree::GeneralFelNode, model::BranchModel; partition_list = 1:length(tree.message), temp_message = deepcopy(tree.message))
    felsenstein_down!(tree, x -> [model], partition_list = partition_list, temp_message = temp_message)
end