begin
    #A mixed partition example
    tree = sim_tree(n = 10)
    nuc_partition = NucleotidePartition(10)
    nuc_partition.state .= 0.25
    internal_message_init!(tree, [nuc_partition, GaussianPartition()])
    Q = HKY85(4.0, 0.25, 0.25, 0.25, 0.25) .* 0.001
    bm_models = [DiagonalizedCTMC(Q), BrownianMotion(0.0, 0.1)]
    sample_down!(tree, bm_models)
    log_likelihood!(tree, bm_models)

    #Testing partition lists
    @test log_likelihood!(tree, bm_models) ≈
          log_likelihood!(tree, bm_models, partition_list = [1]) +
          log_likelihood!(tree, bm_models, partition_list = [2])

    d = marginal_state_dict(tree, bm_models)
    @test typeof.(first(d)[2]) == typeof.(tree.message)
    d = marginal_state_dict(tree, bm_models, partition_list = [1])
    @test typeof.(first(d)[2]) == typeof.(tree.message[[1]])
    d = marginal_state_dict(tree, bm_models, partition_list = [2])
    @test typeof.(first(d)[2]) == typeof.(tree.message[[2]])

    d = cascading_max_state_dict(tree, bm_models)
    @test typeof.(first(d)[2]) == typeof.(tree.message)
    d = cascading_max_state_dict(tree, bm_models, partition_list = [1])
    @test typeof.(first(d)[2]) == typeof.(tree.message[[1]])
    d = cascading_max_state_dict(tree, bm_models, partition_list = [2])
    @test typeof.(first(d)[2]) == typeof.(tree.message[[2]])

    d = endpoint_conditioned_sample_state_dict(tree, bm_models)
    @test typeof.(first(d)[2]) == typeof.(tree.message)
    d = endpoint_conditioned_sample_state_dict(tree, bm_models, partition_list = [1])
    @test typeof.(first(d)[2]) == typeof.(tree.message[[1]])
    d = endpoint_conditioned_sample_state_dict(tree, bm_models, partition_list = [2])
    @test typeof.(first(d)[2]) == typeof.(tree.message[[2]])

    sample_down!(tree, bm_models, partition_list = [1])
    sample_down!(tree, bm_models, partition_list = [2])
    sample_down!(tree, bm_models)
    sample_down!(tree, x -> bm_models, partition_list = [2])
    sample_down!(tree, x -> bm_models)

    log_likelihood!(tree, bm_models, partition_list = [1])
    log_likelihood!(tree, bm_models, partition_list = [2])
    log_likelihood!(tree, bm_models)
    log_likelihood!(tree, x -> bm_models, partition_list = [2])
    log_likelihood!(tree, x -> bm_models)

    felsenstein_down!(tree, bm_models, partition_list = [1])
    felsenstein_down!(tree, bm_models, partition_list = [2])
    felsenstein_down!(tree, bm_models)
    felsenstein_down!(tree, x -> bm_models, partition_list = [2])
    felsenstein_down!(tree, x -> bm_models)

    tol = 1e-5
    unopt_tree_copy = deepcopy(tree)
    branchlength_optim!(tree, bm_models, tol=tol)
    branchlength_sample_under_gss = tree.children[1].branchlength
    branchlength_optim!(unopt_tree_copy, bm_models, tol=tol, bl_optimizer=BrentsMethod())
    branchlength_sample_under_brent = unopt_tree_copy.children[1].branchlength
    @test isapprox(
        MolecularEvolution.unit_inv_transform(branchlength_sample_under_gss), #We're optimizing the untransformed values in the unit domain
        MolecularEvolution.unit_inv_transform(branchlength_sample_under_brent), 
        atol=4*tol) #Brent ensures an accuracy of 3*tol in this case, which is why we choose 4*tol.
    branchlength_optim!(tree, bm_models, partition_list = [1])
    branchlength_optim!(tree, bm_models, partition_list = [2])
    branchlength_optim!(tree, x -> bm_models, partition_list = [2])
    branchlength_optim!(tree, x -> bm_models)

    felsenstein!(tree, bm_models, partition_list = [1])
    felsenstein!(tree, bm_models, partition_list = [2])
    felsenstein!(tree, bm_models)
    felsenstein!(tree, x -> bm_models, partition_list = [2])
    felsenstein!(tree, x -> bm_models)

    nni_optim!(tree, bm_models, partition_list = [1])
    nni_optim!(tree, bm_models, partition_list = [2])
    nni_optim!(tree, bm_models)
    nni_optim!(tree, x -> bm_models, partition_list = [2])
    nni_optim!(tree, x -> bm_models)
end
