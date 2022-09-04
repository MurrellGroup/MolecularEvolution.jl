using MolecularEvolution
using Test

@testset "MolecularEvolution.jl" begin
    @testset "models" begin
        include("test_models.jl")
    end

    @testset "utils" begin
        include("test_utils.jl")
    end

    @testset "optim" begin
        include("test_optim.jl")
    end

    @testset "abstract" begin
        include("test_abstract.jl")
    end

    @testset "felsenstein" begin
        include("test_felsenstein.jl")
    end

    @testset "partition_selection" begin
        include("partition_selection.jl")
    end

end

