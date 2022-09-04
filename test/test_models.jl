begin
    n = GeneralFelNode()
    n.branchlength = 0.1
    Q = [-2.78361 1.22755 0.850393 0.705665; 1.22755 -2.97174 0.96698 0.777209; 0.850393 0.96698 -2.25136 0.433984; 0.705665 0.777209 0.433984 -1.91686]
    m = DiagonalizedCTMC(Q)
    @test sum(abs,P_from_diagonalized_Q(m,n) .- [0.767961 0.0981107 0.0724612 0.061467; 0.0981107 0.754992 0.080452 0.0664456; 0.0724612 0.080452 0.806096 0.0409902; 0.061467 0.0664456 0.0409902 0.831097]) < 0.00001
end

begin
    n = GeneralFelNode()
    n.branchlength = 0.1
    Q = [-2.78361 1.22755 0.850393 0.705665; 1.22755 -2.97174 0.96698 0.777209; 0.850393 0.96698 -2.25136 0.433984; 0.705665 0.777209 0.433984 -1.91686]
    m = DiagonalizedCTMC(Q)
    s = CustomDiscretePartition(4,2)
    d = CustomDiscretePartition(4,2)
    s.state = [0.01 0.02 0.03 0.04; 0.1 0.2 0.3 0.4]'
    backward!(d,s,m,n)
    @test sum(abs,d.state .- [0.0142743 0.142743; 0.0211523 0.211523; 0.0281561 0.281561; 0.0364172 0.364172])<0.00001
end

begin
    n = GeneralFelNode()
    n.branchlength = 1.0
    Q = [-1 1; 0.0 -0.0];
    m = DiagonalizedCTMC(Q)
    s = CustomDiscretePartition(2,3)
    d = CustomDiscretePartition(2,3)
    s.state = [0.1 0.2 0.3; 0.01 0.02 0.03]
    forward!(d,s,m,n)
    @test sum(abs,d.state .- [0.0367879 0.0735759 0.110364; 0.0732121 0.146424 0.219636]) < 0.000001
end

