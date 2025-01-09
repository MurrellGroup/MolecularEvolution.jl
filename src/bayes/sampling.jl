"""
    metropolis_sample(
        update!::Function,
        initial_tree::FelNode,
        models::Vector{<:BranchModel},
        num_of_samples;
        burn_in = 1000,
        sample_interval = 10,
        collect_LLs = false,
        midpoint_rooting = false,
        ladderize = false,
    )

Samples tree topologies from a posterior distribution using a custom `update!` function.

# Arguments
- `update!`: A function that takes (tree::FelNode, models::Vector{<:BranchModel}) and updates `tree`. `update!` takes (tree::FelNode, models::Vector{<:BranchModel}) and updates `tree`. One call to `update!` corresponds to one iteration of the Metropolis algorithm.
- `initial_tree`: An initial tree topology with the leaves populated with data, for the likelihood calculation.
- `models`: A list of branch models.
- `num_of_samples`: The number of tree samples drawn from the posterior.
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
    update!::Function,
    initial_tree::FelNode,
    models::Vector{<:BranchModel},
    num_of_samples;
    burn_in = 1000,
    sample_interval = 10,
    collect_LLs = false,
    midpoint_rooting = false,
    ladderize = false,
)

    # The prior over the (log) of the branchlengths should be specified in bl_sampler. 
    # Furthermore, a non-informative/uniform prior is assumed over the tree topolgies (excluding the branchlengths).

    sample_LLs = []
    samples = FelNode[]
    tree = deepcopy(initial_tree)
    iterations = burn_in + num_of_samples * sample_interval

    for i = 1:iterations
        # Updates the tree topolgy and branchlengths.
        update!(tree, models)

        if (i - burn_in) % sample_interval == 0 && i > burn_in

            push!(samples, copy_tree(tree, true))

            if collect_LLs
                push!(sample_LLs, log_likelihood!(tree, models))
            end

        end

    end

    if midpoint_rooting
        for (i, sample) in enumerate(samples)
            node, len = midpoint(sample)
            samples[i] = reroot!(node, dist_above_child = len)
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

"""
    metropolis_sample(
        initial_tree::FelNode,
        models::Vector{<:BranchModel},
        num_of_samples;
        bl_sampler::UnivariateSampler = BranchlengthSampler(Normal(0,2), Normal(-1,1))
        burn_in=1000, 
        sample_interval=10,
        collect_LLs = false,
        midpoint_rooting=false,
    )

A convenience method. One step of the Metropolis algorithm is performed by calling [`nni_update!`](@ref) with `softmax_sampler` and [`branchlength_update!`](@ref) with `bl_sampler`.

# Additional Arguments
- `bl_sampler`: Sampler used to drawn branchlengths from the posterior.
"""
function metropolis_sample(
    args...;
    bl_sampler::UnivariateSampler = BranchlengthSampler(Normal(0, 2), Normal(-1, 1)),
    kwargs...,
)
    metropolis_sample(args...; kwargs...) do tree, models
        nni_update!(softmax_sampler, tree, x -> models)
        branchlength_update!(bl_sampler, tree, x -> models)
    end
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

    sort!(leaflist, by = x -> x.name)

    order = [node_dic[leaf] for leaf in leaflist]

    return distmat[order, order]
end
