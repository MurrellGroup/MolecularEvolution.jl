
export metropolis_sample
"""
    function metropolis_sample(
        initial_tree::FelNode,
        models::Vector{<:BranchModel},
        num_of_samples;
        burn_in=1000, 
        sample_interval=10,
        collect_LLs = false,
        midpoint_rooting=false,
    )

Samples tree topologies from a posterior distribution.

# Arguments
- `initial_tree`: An initial topology with (important!) the leaves populated with data, for the likelihood calculation.
- `models`: A list of branch models.
- `num_of_samples`: The number of tree samples drawn from the posterior.
- `burn_in`: The number of samples discarded at the start of the Markov Chain.
- `sample_interval`: The distance between samples in the underlying Markov Chain (to reduce sample correlation).
- `collect_LLs`: Specifies if the function should return the log-likelihoods of the trees.
- `midpoint_rooting`: Specifies whether the drawn samples should be midpoint rerooted (Important! Should only be used for time-reversible branch models starting in equilibrium).

# Returns
- `samples`: The trees drawn from the posterior.
- `sample_LLs`: The associated log-likelihoods of the tree (optional).
"""
function metropolis_sample(
    initial_tree::FelNode,
    models::Vector{<:BranchModel},
    num_of_samples;
    burn_in=1000, 
    sample_interval=10,
    collect_LLs = false,
    midpoint_rooting=false,
    ladderize = false,
    )

    sample_LLs = []
    samples = FelNode[]
    tree = deepcopy(initial_tree)
    iterations = burn_in + num_of_samples * sample_interval

    bl_modifier = BranchlengthPerturbation(2.0)

    #old_softmax_sampler = x -> (sample = rand(Categorical(softmax(x))); changed = sample != 1; (changed, sample))
    softmax_sampler = x -> rand(Categorical(softmax(x)))
    for i=1:iterations

            # Updates the tree topolgy and branchlengths using Gibbs sampling.
            nni_optim!(tree, models, nni_config_sampler = softmax_sampler)
            branchlength_optim!(tree, models, modifier=bl_modifier)   

            if (i-burn_in) % sample_interval == 0 && i > burn_in

                push!(samples, copy_tree(tree, true))

                if collect_LLs
                    push!(sample_LLs, log_likelihood!(tree, models))
                end
                
            end

            #REMOVE BEFORE PR
            if i % 1000 == 0 || i == iterations 
                println(floor(i/iterations * 100))
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

function branchlength_metropolis(LL, modifier, curr_value)
    # The prior distribution for the variable log(branchlength). A small perturbation of +1e-12 is added to enhance numerical stability near zero.
    log_prior(x) = logpdf(Normal(-1,1),log(x + 1e-12))
    # Adding additive normal symmetrical noise in the log(branchlength) domain to ensure the proposal function is symmetric.
    noise = modifier.sigma*rand(Normal(0,1))
    proposal = exp(log(curr_value)+noise)
    # The standard Metropolis acceptance criterion.
    if rand() <= exp(LL(proposal)+log_prior(proposal)-LL(curr_value)-log_prior(curr_value))
        modifier.acc_ratio[1] = modifier.acc_ratio[1] + 1
        return proposal
    else
        modifier.acc_ratio[2] = modifier.acc_ratio[2] + 1 
        return curr_value
    end
end

export collect_leaf_dists
"""
    collect_leaf_dists(trees::Vector{<:AbstractTreeNode})
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


