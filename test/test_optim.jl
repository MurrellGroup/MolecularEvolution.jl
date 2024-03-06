begin
    f(x) = -(x - 2)^2
    m = golden_section_maximize(f, 1, 5, identity, 1e-20)
    @test m == 2.0
    m = brents_method_minimize(x -> -f(x), 1, 5, identity, 1e-20)
    @test m == 2.0
end

begin
    @test unit_transform(0.5) == 1.0
end
