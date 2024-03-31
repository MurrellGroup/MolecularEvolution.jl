begin
    n = FelNode()
    n.branchlength = 0.1
    Q = [
        -2.78361 1.22755 0.850393 0.705665
        1.22755 -2.97174 0.96698 0.777209
        0.850393 0.96698 -2.25136 0.433984
        0.705665 0.777209 0.433984 -1.91686
    ]
    m = DiagonalizedCTMC(Q)
    @test sum(
        abs,
        P_from_diagonalized_Q(m, n) .- [
            0.767961 0.0981107 0.0724612 0.061467
            0.0981107 0.754992 0.080452 0.0664456
            0.0724612 0.080452 0.806096 0.0409902
            0.061467 0.0664456 0.0409902 0.831097
        ],
    ) < 0.00001
end

begin
    n = FelNode()
    n.branchlength = 0.1
    Q = [
        -2.78361 1.22755 0.850393 0.705665
        1.22755 -2.97174 0.96698 0.777209
        0.850393 0.96698 -2.25136 0.433984
        0.705665 0.777209 0.433984 -1.91686
    ]
    m = DiagonalizedCTMC(Q)
    s = CustomDiscretePartition(4, 2)
    d = CustomDiscretePartition(4, 2)
    s.state = [0.01 0.02 0.03 0.04; 0.1 0.2 0.3 0.4]'
    backward!(d, s, m, n)
    @test sum(
        abs,
        d.state .-
        [0.0142743 0.142743; 0.0211523 0.211523; 0.0281561 0.281561; 0.0364172 0.364172],
    ) < 0.00001
end

begin
    n = FelNode()
    n.branchlength = 1.0
    Q = [-1 1; 0.0 -0.0]
    m = DiagonalizedCTMC(Q)
    s = CustomDiscretePartition(2, 3)
    d = CustomDiscretePartition(2, 3)
    s.state = [0.1 0.2 0.3; 0.01 0.02 0.03]
    forward!(d, s, m, n)
    @test sum(abs, d.state .- [0.0367879 0.0735759 0.110364; 0.0732121 0.146424 0.219636]) <
          0.000001
end

begin #PiQ
    n = FelNode()
    n.branchlength = 1.0
    dest = NucleotidePartition(10)
    dest2 = NucleotidePartition(10)
    src = NucleotidePartition(10)
    src.state .= rand(4,10)

    m = PiQ([0.1,0.2,0.3,0.4])

    #Now construct a Q matrix that matches this model, where each column is equal to the PiQ freqs, except the diagonals which must be the negative sum of the rest of the rows
    Q = [m.pi[j] for i in 1:4, j in 1:4]
    for i in 1:4
        Q[i,i] = -sum(Q[i,j] for j in 1:4 if j != i)
    end
    m2 = GeneralCTMC(Q)
    
    forward!(dest,src,m,n)
    forward!(dest2,src,m2,n)

    @test dest2.state ≈ dest.state

    backward!(dest,src,m,n)
    backward!(dest2,src,m2,n)

    @test dest2.state ≈ dest.state

    n.branchlength = 1000.0
    forward!(dest,src,m,n)
    forward!(dest2,src,m2,n)

    @test dest2.state ≈ dest.state

    n.branchlength = 1000.0
    backward!(dest,src,m,n)
    backward!(dest2,src,m2,n)

    @test dest2.state ≈ dest.state
end

begin #LazyPartition
    """
    n = FelNode()
    n.branchlength = 1.0
    #Create some CodonPartitions
    sites = 10
    codon_dest = CodonPartition(sites)
    codon_dest.state .= rand(length(MolecularEvolution.universal_code.sense_codons))
    codon_src = CodonPartition(sites)
    codon_src.state .= rand(length(MolecularEvolution.universal_code.sense_codons))

    #Wrap copies of them in LazyPartitions
    memoryblocks = Vector{CodonPartition}()
    lazy_dest = LazyPartition{CodonPartition}(deepcopy(codon_dest))
    lazy_src = LazyPartition{CodonPartition}(deepcopy(codon_src))
    lazy_dest.memoryblocks, lazy_src.memoryblocks = memoryblocks

    #Check that the state is the same after a combine!
    combine!(codon_dest, codon_src)
    combine!(lazy_dest, lazy_src)

    @test codon_dest.state ≈ lazy_dest.partition.state
    @test isnothing(lazy_src.partition)
    
    #Switch order such that lazy_dest.partition is nothing
    codon_dest, codon_src = codon_src, codon_dest
    lazy_dest, lazy_src = lazy_src, lazy_dest

    #Check that the state is the same after a backward!
    model = GeneralCTMC(MolecularEvolution.MG94_F3x4(1, 1, reversibleQ([1.0,2.0,1.0,1.0,1.6,0.5],[0.2,0.3,0.3,0.2]), [0.2 0.3 0.3 0.2; 0.25 0.25 0.25 0.25; 0.3 0.2 0.2 0.3]))
    backward!(codon_dest, codon_src, model, n)
    backward!(lazy_dest, lazy_src, model, n)

    @test codon_dest.state ≈ lazy_dest.partition.state
    """

    unopt_tree = sim_tree(n=15)
    opt_tree = deepcopy(unopt_tree)

    maximum_active_partitions = MolecularEvolution.lazysort!(opt_tree)
    #@show maximum_active_partitions, treedepth(unopt_tree), length(getnodelist(unopt_tree))

    @test maximum_active_partitions <= treedepth(unopt_tree)
    @test maximum_active_partitions <= floor(log2(length(getnodelist(unopt_tree)))) + 1
    
    """
    p1 = plot(get_phylo_tree(unopt_tree), size=(400, 400))
    p2 = plot(get_phylo_tree(opt_tree), size=(400,400))
    plot(p1, p2, layout=(2,1))
    """
    #TODO test this on a tree! Both speed, allocs memory usage.
end