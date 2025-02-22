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

    #TODO When we use BrentsMethodOpt, check if we gain a speed-up and that we're not catastrophically wrong
    branchlength_optim!(tree, bm_models, partition_list = [1])
    branchlength_optim!(tree, bm_models, partition_list = [2])
    branchlength_optim!(tree, bm_models)
    branchlength_optim!(tree, bm_models, bl_optimizer=BrentsMethodOpt())
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

    felsenstein_roundtrip!(tree, bm_models, partition_list = [1])
    felsenstein_roundtrip!(tree, bm_models, partition_list = [2])
    felsenstein_roundtrip!(tree, bm_models)
    felsenstein_roundtrip!(tree, x -> bm_models, partition_list = [2])
    felsenstein_roundtrip!(tree, x -> bm_models)

    tree = root_optim!(tree, bm_models, partition_list = [1])
    tree = root_optim!(tree, bm_models, partition_list = [2])
    tree = root_optim!(tree, bm_models)
    tree = root_optim!(tree, x -> bm_models, partition_list = [2])
    tree = root_optim!(tree, x -> bm_models)

    updater = MaxLikUpdate(root=1)
    tree, bm_models = updater(tree, bm_models, partition_list = [1])
    updater = BayesUpdate(root=0, models=1)
    tree, bm_models = updater(tree, bm_models, partition_list = [1])
end
