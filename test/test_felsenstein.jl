begin
    newt = gettreefromnewick("((A:0.023,B:0.234):0.124,(C:0.123,D:0.234):0.324);", FelNode)
    Q = HKY85(5.0, 0.1, 0.2, 0.3, 0.4)
    m1 = DiagonalizedCTMC(Q)
    function models(node::FelNode)
        return [m1]
    end
    empty_mess = [NucleotidePartition(2)]
    newt.parent_message = equilibrium_message([m1], empty_mess)
    internal_message_init!(newt, empty_mess)
    for n in getnodelist(newt)
        n.message = deepcopy(empty_mess)
    end
    for n in getnodelist(newt)
        n.message[1].state = [1.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 1.0]
    end
    n = getleaflist(newt)[1]
    n.message[1].state = [0.0 1.0; 0.0 0.0; 0.0 0.0; 1.0 0.0]
    null_mess = deepcopy(newt.message)
    newt.parent_message = equilibrium_message([m1], empty_mess)
    # TODO: The following code results in a method error.
    #felsenstein_up!(null_mess,newt,models,1:length(null_mess))
    #@test log_likelihood(newt) == -18.05871954142185 #Needs checking!!
    #@test sum(abs2,newt.message[1].state .- [0.985326 0.42005; 0.00055197 0.0085463; 0.0010581 0.00453591; 0.0130636 0.566867]) < 0.0000001
end

begin
    s = GappyNucleotidePartition(10)
    s.state = rand(size(s.state)...)
    sample_partition!(s)
    @test sum([sum(s.state[:, i] .^ 2) for i = 1:10]) == 10.0
end

begin
    d = CustomDiscretePartition(2, 2)
    s_arr = [CustomDiscretePartition(2, 2), CustomDiscretePartition(2, 2)]
    s_arr[1].state, s_arr[2].state = ([1 2; 3 4], [10 20; 30 40])
    combine!(d, s_arr, true)
    @test d.state == [0.1 0.2; 0.9 0.8]
end
