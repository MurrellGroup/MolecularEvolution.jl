let #Rec. vs. Iter.
    partlists = [1:1, 2:2, 1:2]
    for i = 1:500
        partlist = partlists[mod(i, 3)+1]
        tree = sim_tree(n=rand(20:500))
        nuc_partition = NucleotidePartition(10)
        nuc_partition.state .= 0.25
        internal_message_init!(tree, [nuc_partition, GaussianPartition()])
        Q = HKY85(4.0, 0.25, 0.25, 0.25, 0.25) .* 0.001
        bm_models = [DiagonalizedCTMC(Q), BrownianMotion(0.0, 1.0)]
        sample_down!(tree, bm_models)
        for n in getnodelist(tree)
            n.branchlength *= (rand()+0.5)
        end
        felsenstein!(tree, bm_models, partition_list = partlist)
        tree_copy = deepcopy(tree)
        branchlength_optim!(tree, bm_models, partition_list = partlist)
        MolecularEvolution.branchlength_optim_iter!(tree_copy, x -> bm_models, partition_list = partlist)

        branchlengths = [node.branchlength for node in getnodelist(tree)]
        branchlengths_copy = [node.branchlength for node in getnodelist(tree_copy)]
        
        @assert isapprox(branchlengths, branchlengths_copy)
        branchlengths != branchlengths_copy && @warn ""
    end
end