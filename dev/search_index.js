var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = MolecularEvolution","category":"page"},{"location":"#MolecularEvolution","page":"Home","title":"MolecularEvolution","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for MolecularEvolution.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [MolecularEvolution]","category":"page"},{"location":"#Base.:==-Union{Tuple{T}, Tuple{T, T}} where T<:AbstractTreeNode","page":"Home","title":"Base.:==","text":"==(t1, t2)\nDefaults to pointer equality\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution._mapreduce-Union{Tuple{T}, Tuple{AbstractTreeNode, T, Any, Any}} where T<:Function","page":"Home","title":"MolecularEvolution._mapreduce","text":"Internal function. Helper for bfsmapreduce and dfsmapreduce\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.bfs_mapreduce-Union{Tuple{T}, Tuple{AbstractTreeNode, T, Any}} where T<:Function","page":"Home","title":"MolecularEvolution.bfs_mapreduce","text":"Performs a BFS map-reduce over the tree, starting at a given node For each node, mapreduce is called as:    mapreduce(currnode::GeneralFelNode, prevnode::GeneralFelNode, aggregator) where prev_node is the previous node visited on the path from the start node to the current node It is expected to update the aggregator, and not return anything.\n\nNot exactly conventional map-reduce, as map-reduce calls may rely on state in the aggregator added by map-reduce calls on other nodes visited earlier.\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.branchlength_optim!-Tuple{GeneralFelNode, Any}","page":"Home","title":"MolecularEvolution.branchlength_optim!","text":"branchlength_optim!(tree::GeneralFelNode, models; partition_list = nothing, tol = 1e-5)\n\nUses golden section search to optimize all branches recursively, maintaining the integrity of the messages. Requires felsenstein!() to have been run first.\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.deepequals-Union{Tuple{T}, Tuple{T, T}} where T<:AbstractTreeNode","page":"Home","title":"MolecularEvolution.deepequals","text":"deepequals(t1, t2)\n\nChecks whether two trees are equal by recursively calling this on all fields, except :parent, in order to prevent cycles. In order to ensure that the :parent field is not hiding something different on both trees, ensure that each is consistent first (see: istreeconsistent).\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.dfs_mapreduce-Union{Tuple{T}, Tuple{AbstractTreeNode, T, Any}} where T<:Function","page":"Home","title":"MolecularEvolution.dfs_mapreduce","text":"Performs a DFS map-reduce over the tree, starting at a given node See bfs_mapreduce for more details.\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.expected_subs_per_site-Tuple{Any, Any}","page":"Home","title":"MolecularEvolution.expected_subs_per_site","text":"expected_subs_per_site(Q,mu)\n\nTakes a rate matrix Q and an equilibrium frequency vector, and calculates the expected number of substitutions per site.\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.felsenstein!-Tuple{GeneralFelNode, Any}","page":"Home","title":"MolecularEvolution.felsenstein!","text":"felsenstein!(node::GeneralFelNode, models; partition_list = nothing)\n\nShould usually be called on the root of the tree. Propagates Felsenstein pass up from the tips to the root. models must be a function that takes a node, and returns a Vector{BranchModel}. This lets you control which models gets used for which branch. partition_list (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over.\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.felsenstein!-Tuple{GeneralFelNode, BranchModel}","page":"Home","title":"MolecularEvolution.felsenstein!","text":"felsenstein!(node::GeneralFelNode, model::BranchModel; partition_list = nothing)\n\nShould usually be called on the root of the tree. Propagates Felsenstein pass up from the tips to the root. models can be a single model, or a vector of models, or a function that associates a node with a model. partition_list (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over.\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.felsenstein!-Tuple{GeneralFelNode, Vector{<:BranchModel}}","page":"Home","title":"MolecularEvolution.felsenstein!","text":"felsenstein!(node::GeneralFelNode, models::Vector{<:BranchModel}; partition_list = nothing)\n\nShould usually be called on the root of the tree. Propagates Felsenstein pass up from the tips to the root. models can be a single model, or a vector of models, or a function that associates a node with a model. partition_list (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over.\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.felsenstein_down!-Tuple{GeneralFelNode, Any}","page":"Home","title":"MolecularEvolution.felsenstein_down!","text":"felsenstein_down!(node::GeneralFelNode, models; partition_list = 1:length(tree.message), temp_message = deepcopy(tree.message))\n\nShould usually be called on the root of the tree. Propagates Felsenstein pass down from the root to the tips. felsenstein!() should usually be called first. models must be a function that takes a node, and returns a Vector{BranchModel}. This lets you control which models gets used for which branch. partition_list (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over.\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.felsenstein_down!-Tuple{GeneralFelNode, BranchModel}","page":"Home","title":"MolecularEvolution.felsenstein_down!","text":"felsenstein_down!(node::GeneralFelNode, models; partition_list = 1:length(tree.message), temp_message = deepcopy(tree.message))\n\nShould usually be called on the root of the tree. Propagates Felsenstein pass down from the root to the tips. felsenstein!() should usually be called first. models can be a single model, or a vector of models, or a function that associates a node with a model. partition_list (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over.\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.felsenstein_down!-Tuple{GeneralFelNode, Vector{<:BranchModel}}","page":"Home","title":"MolecularEvolution.felsenstein_down!","text":"felsenstein_down!(node::GeneralFelNode, models; partition_list = 1:length(tree.message), temp_message = deepcopy(tree.message))\n\nShould usually be called on the root of the tree. Propagates Felsenstein pass down from the root to the tips. felsenstein!() should usually be called first. models can be a single model, or a vector of models, or a function that associates a node with a model. partition_list (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over.\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.golden_section_maximize-Tuple{Any, Real, Real, Any, Real}","page":"Home","title":"MolecularEvolution.golden_section_maximize","text":"Golden section search.\n\nGiven a function f with a single local minimum in the interval [a,b], gss returns a subset interval [c,d] that contains the minimum with d-c <= tol.\n\nExamples\n\njulia> f(x) = -(x-2)^2\nf (generic function with 1 method)\n\njulia> m = golden_section_maximize(f, 1, 5, identity, 1e-10)\n2.0000000000051843\n\nFrom: https://en.wikipedia.org/wiki/Golden-section_search\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.internal_message_init!-Tuple{GeneralFelNode, Partition}","page":"Home","title":"MolecularEvolution.internal_message_init!","text":"internal_message_init!(tree::GeneralFelNode, partition::Partition)\n\nInitializes the message template for each node in the tree, as an array of the partition.\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.internal_message_init!-Tuple{GeneralFelNode, Vector{<:Partition}}","page":"Home","title":"MolecularEvolution.internal_message_init!","text":"internal_message_init!(tree::GeneralFelNode, empty_message::Vector{<:Partition})\n\nInitializes the message template for each node in the tree, allocating space for each partition.\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.istreeconsistent-Tuple{T} where T<:AbstractTreeNode","page":"Home","title":"MolecularEvolution.istreeconsistent","text":"istreeconsistent(root)\n\nChecks whether the :parent field is set to be consistent with the :child field for all nodes in the subtree. \n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.linear_scale-NTuple{5, Any}","page":"Home","title":"MolecularEvolution.linear_scale","text":"linear_scale(val,in_min,in_max,out_min,out_max)\n\nLinearly maps val which lives in [inmin,inmax] to a value in [outmin,outmax]\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.log_likelihood!-Tuple{GeneralFelNode, Any}","page":"Home","title":"MolecularEvolution.log_likelihood!","text":"log_likelihood!(tree::GeneralFelNode, models; partition_list = nothing)\n\nFirst re-computes the upward felsenstein pass, and then computes the log likelihood of this tree.\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.log_likelihood-Tuple{GeneralFelNode, Any}","page":"Home","title":"MolecularEvolution.log_likelihood","text":"log_likelihood(tree::GeneralFelNode; partition_list = nothing)\n\nComputed the log likelihood of this tree. Requires felsenstein!() to have been run.\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.longest_path-Tuple{GeneralFelNode}","page":"Home","title":"MolecularEvolution.longest_path","text":"Returns the longest path in a tree For convenience, this is returned as two lists of form:     [leafnode, parentnode, .... root] Where the leaf_node nodes are selected to be the furthest away\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.midpoint-Tuple{GeneralFelNode}","page":"Home","title":"MolecularEvolution.midpoint","text":"Returns a midpoint as a node and a distance above it where the midpoint is\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.name2node_dict-Tuple{AbstractTreeNode}","page":"Home","title":"MolecularEvolution.name2node_dict","text":"name2node_dict(root)\n\nReturns a dictionary of leaf nodes, indexed by node.name. Can be used to associate sequences with leaf nodes.\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.node_distances-Tuple{AbstractTreeNode}","page":"Home","title":"MolecularEvolution.node_distances","text":"Compute the distance to all other nodes from a given node\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.parent_list-Tuple{GeneralFelNode}","page":"Home","title":"MolecularEvolution.parent_list","text":"Provides a list of parent nodes nodes from this node up to the root node\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.quadratic_CI-Tuple{Function, Vector, Int64}","page":"Home","title":"MolecularEvolution.quadratic_CI","text":"quadratic_CI(f::Function,opt_params::Vector, param_ind::Int; rate_conf_level = 0.99, nudge_amount = 0.01)\n\nTakes a NEGATIVE log likelihood function (compatible with Optim.jl), a vector of maximizing parameters, an a parameter index. Returns the quadratic confidence interval.\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.quadratic_CI-Tuple{Vector, Vector}","page":"Home","title":"MolecularEvolution.quadratic_CI","text":"quadratic_CI(xvec,yvec; rate_conf_level = 0.99)\n\nTakes xvec, a vector of parameter values, and yvec, a vector of log likelihood evaluations (note: NOT the negative LLs you) might use with Optim.jl. Returns the confidence intervals computed by a quadratic approximation to the LL.\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.root2tip_distances-Tuple{AbstractTreeNode}","page":"Home","title":"MolecularEvolution.root2tip_distances","text":"root2tips(root::AbstractTreeNode)\n\nReturns a vector of root-to-tip distances, and a node-to-index dictionary. Be aware that this dictionary will break when any of the node content (ie. anything on the tree) changes.\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.sample_down!","page":"Home","title":"MolecularEvolution.sample_down!","text":"sampledown!(node::GeneralFelNode,models::Vector{<:BranchModel},partitionlist)\n\nGenerates samples under the model, matching the partition structure.\n\n\n\n\n\n","category":"function"},{"location":"#MolecularEvolution.sample_down!-Tuple{GeneralFelNode, Any, Any}","page":"Home","title":"MolecularEvolution.sample_down!","text":"sampledown!(node::GeneralFelNode,models,partitionlist)\n\nGenerates samples under the model. The root.parent_message is taken as the starting distribution, and node.message contains the sampled messages.\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.sample_down!-Tuple{GeneralFelNode, Any}","page":"Home","title":"MolecularEvolution.sample_down!","text":"sample_down!(node::GeneralFelNode,models)\n\nGenerates samples under the model, matching the partition structure.\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.sample_down!-Tuple{GeneralFelNode, BranchModel}","page":"Home","title":"MolecularEvolution.sample_down!","text":"sampledown!(node::GeneralFelNode,model::BranchModel,partitionlist)\n\nGenerates samples under the model, matching the partition structure.\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.shortest_path_between_nodes-Tuple{GeneralFelNode, GeneralFelNode}","page":"Home","title":"MolecularEvolution.shortest_path_between_nodes","text":"Shortest path between nodes, returned as two lists, each starting with one of the two nodes,  and ending with the common ancestor\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.sibling_inds-Tuple{AbstractTreeNode}","page":"Home","title":"MolecularEvolution.sibling_inds","text":"sibling_inds(node)\n\nReturns logical indices of the siblings in the parent's child's vector.\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.siblings-Tuple{AbstractTreeNode}","page":"Home","title":"MolecularEvolution.siblings","text":"siblings(node)\n\nReturns a vector of siblings of node.\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.sim_tree-Tuple{Int64, Any, Any}","page":"Home","title":"MolecularEvolution.sim_tree","text":"sim_tree(add_limit::Int,Ne_func,sample_rate_func; nstart = 1, time = 0.0, mutation_rate = 1.0, T = Float64)\n\nSimulates a tree of type GeneralFelNode{T}. Allows an effective population size function (Nefunc), as well as a sample rate function (samplerate_func), which can also just be constants.\n\nNefunc(t) = (sin(t/10)+1)*100.0 + 10.0 root = simtree(600,Nefunc,1.0) simpletree_draw(ladderize(root))\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.sim_tree-Tuple{}","page":"Home","title":"MolecularEvolution.sim_tree","text":"sim_tree(;n = 10)\n\nSimulates tree with constant population size.\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.simple_tree_opt!-Tuple{Any, Any}","page":"Home","title":"MolecularEvolution.simple_tree_opt!","text":"simple_tree_opt!(newt,models; tol = 10^-4, verbose = 0, topology = true)\n\nTakes a tree and a model function, and optimizes branch lengths and, optionally, topology.\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.tree2distances-Tuple{AbstractTreeNode}","page":"Home","title":"MolecularEvolution.tree2distances","text":"tree2distances(root::AbstractTreeNode)\n\nReturns a distance matrix for all pairs of leaf nodes, and a node-to-index dictionary. Be aware that this dictionary will break when any of the node content (ie. anything on the tree) changes.\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.tree2shared_branch_lengths-Tuple{AbstractTreeNode}","page":"Home","title":"MolecularEvolution.tree2shared_branch_lengths","text":"tree2distances(root::AbstractTreeNode)\n\nReturns a distance matrix for all pairs of leaf nodes, and a node-to-index dictionary. Be aware that this dictionary will break when any of the node content (ie. anything on the tree) changes.\n\n\n\n\n\n","category":"method"},{"location":"#MolecularEvolution.tree_from_file-Tuple{Any}","page":"Home","title":"MolecularEvolution.tree_from_file","text":"tree_from_file(treefile)\n\nReads in a tree from a file, of type GeneralFelNode\n\n\n\n\n\n","category":"method"}]
}
