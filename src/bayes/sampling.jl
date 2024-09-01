"""
    function metropolis_sample(
        initial_tree::FelNode,
        models::Vector{<:BranchModel},
        num_of_samples;
        bl_modifier::UnivariateSampler = BranchlengthSampler(Normal(0,2), Normal(-1,1))
        burn_in=1000, 
        sample_interval=10,
        collect_LLs = false,
        midpoint_rooting=false,
    )

Samples tree topologies from a posterior distribution. 

# Arguments
- `initial_tree`: An initial tree topology with the leaves populated with data, for the likelihood calculation.
- `models`: A list of branch models.
- `num_of_samples`: The number of tree samples drawn from the posterior.
- `bl_sampler`: Sampler used to drawn branchlengths from the posterior. 
- `burn_in`: The number of samples discarded at the start of the Markov Chain.
- `sample_interval`: The distance between samples in the underlying Markov Chain (to reduce sample correlation).
- `collect_LLs`: Specifies if the function should return the log-likelihoods of the trees.
- `midpoint_rooting`: Specifies whether the drawn samples should be midpoint rerooted (Important! Should only be used for time-reversible branch models starting in equilibrium).

!!! note
    The leaves of the initial tree should be populated with data and felsenstein! should be called on the initial tree before calling this function.

# Returns
- `samples`: The trees drawn from the posterior. Returns shallow tree copies, which needs to be repopulated before running felsenstein! etc. 
- `sample_LLs`: The associated log-likelihoods of the tree (optional).
"""
function metropolis_sample(
    initial_tree::FelNode,
    models::Vector{<:BranchModel},
    num_of_samples;
    bl_sampler::UnivariateSampler = BranchlengthSampler(Normal(0,2), Normal(-1,1)),
    burn_in=1000, 
    sample_interval=10,
    collect_LLs = false,
    midpoint_rooting=false,
    ladderize = false,
    )

    # The prior over the (log) of the branchlengths should be specified in bl_sampler. 
    # Furthermore, a non-informative/uniform prior is assumed over the tree topolgies (excluding the branchlengths).

    sample_LLs = []
    samples = FelNode[]
    tree = deepcopy(initial_tree)
    iterations = burn_in + num_of_samples * sample_interval
    
    softmax_sampler = x -> rand(Categorical(softmax(x)))
    for i=1:iterations
        
            # Updates the tree topolgy and branchlengths.
            nni_optim!(tree, x -> models, selection_rule = softmax_sampler)
            branchlength_optim!(tree, x -> models, bl_modifier = bl_sampler)   

            if (i-burn_in) % sample_interval == 0 && i > burn_in

                push!(samples, copy_tree(tree, true))

                if collect_LLs
                    push!(sample_LLs, log_likelihood!(tree, models))
                end
                
            end

    end

    if midpoint_rooting
        for (i,sample) in enumerate(samples)
            node, len = midpoint(sample)
            samples[i] = reroot!(node, dist_above_child=len)
        end
    end

    if ladderize
        for sample in samples
            ladderize!(sample)
        end
    end

    if collect_LLs
        return samples, sample_LLs
    end

    return samples
end

# Below are some functions that help to assess the mixing by looking at the distance between leaf nodes.

"""
    collect_leaf_dists(trees::Vector{<:AbstractTreeNode})

    Returns a list of distance matrices containing the distance between the leaf nodes, which can be used to assess mixing.
"""
function collect_leaf_dists(trees::Vector{<:AbstractTreeNode})
    distmats = []
    for tree in trees
        push!(distmats, leaf_distmat(tree))
    end
    return distmats
end

"""
    leaf_distmat(tree)

 Returns a matrix of the distances between the leaf nodes where the index on the columns and rows are sorted by the leaf names.
"""
function leaf_distmat(tree)
    
    distmat, node_dic = MolecularEvolution.tree2distances(tree)
    
    leaflist = getleaflist(tree)
    
    sort!(leaflist, by = x-> x.name)
    
    order = [node_dic[leaf] for leaf in leaflist]
    
    return distmat[order, order]
end


