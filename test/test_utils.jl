begin
    testv = [0.1,0.2,0.3,0.4]
    @test sum(abs2,one_hot_sample(testv)) == 1.0
    @test length(one_hot_sample(testv)) == 4
    testv = Vector{Float32}([0.1,0.2,0.3,0.4])
    @test eltype(one_hot_sample(testv)) == eltype(testv)
end

begin
    vec = [-3.0,-7.0]
    scaled_prob_domain(vec) == [1.0, 0.0183156]
end

begin
    m = [0.642979 0.903402 0.757445; 0.384552 0.553417 0.788832; 0.808534 0.811944 0.907218; 0.104483 0.78629 0.427887]
    v = [0.596797, 0.0649148, 0.935419]
    scale_cols_by_vec!(m,v)
    sum(abs2,[0.229008 0.00380687 0.662771; 0.136965 0.00233206 0.690234; 0.287973 0.00342147 0.793823; 0.0372134 0.00331337 0.374405])<0.000001
end
