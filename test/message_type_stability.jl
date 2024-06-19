begin
    # Single partition example
    tree = sim_tree(n = 10)
    GAA_partition = GappyAminoAcidPartition(5)
    AA_freqs = [1 / GAA_partition.states for _ = 1:GAA_partition.states]
    GAA_partition.state .= AA_freqs
    internal_message_init!(tree, GAA_partition)
    Q = gappy_Q_from_symmetric_rate_matrix(WAGmatrix, 1.0, AA_freqs)
    model = DiagonalizedCTMC(Q)
    sample_down!(tree, model)
    felsenstein!(tree, model)

    # These would previously break since Vector{GappyAminoAcidPartition} is not <: Vector{Partition}, for example.
    branchlength_optim!(tree, model)
    marginal_state_dict(tree, model)
    nni_optim!(tree, model)
end