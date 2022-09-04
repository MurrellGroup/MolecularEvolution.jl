begin
    newt = gettreefromnewick("((A:0.023,B:0.234):0.124,(C:0.123,D:0.234):0.324);", GeneralFelNode);
    @test length(getnodelist(newt)) == 7
    @test sum([node.branchlength for node in getnodelist(newt)]) == 1.062
end

begin
    newt = gettreefromnewick("((A:0.023,B:0.234):0.124,(C:0.123,D:0.234):0.324);", GeneralFelNode);
    midpoint_node, distance_above = midpoint(newt)
    new_tree = reroot!(midpoint_node; dist_above_child = distance_above)
    longest_path_1, longest_path_2 = longest_path(new_tree)
    @test sum(x->x.branchlength, longest_path_1) ≈ sum(x->x.branchlength, longest_path_2) 
end
