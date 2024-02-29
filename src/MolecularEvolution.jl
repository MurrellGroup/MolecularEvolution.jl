module MolecularEvolution
#This is rougly divided up into five parts:
#1. "core": which contains the tree code, the message passing algorithms, and the tree topology/branch length optimization code
#2. "models": containing type defs for types of states, and models that act on them
#3. "optimization": containing code for parameter optimization - underdeveloped
#4. "bayes": containing code for bayesian inference - underdeveloped
#5. "plots": containing code for drawing trees

#Whenever possible, we'll keep dependencies to a minimum using Requires.jl
using Requires

using Distributions
using LinearAlgebra
using StatsBase

#Internal constants
const eps = 1e-15
#const e = MathConstants.e


abstract type Partition end
abstract type MultiSitePartition <: Partition end #Requires partition.scaling to be an array of log likelihoods, one per site
abstract type DiscretePartition <: MultiSitePartition end #For convenience, all our DiscretePartitions are MultiSite
abstract type ContinuousPartition <: Partition end

abstract type BranchModel end
abstract type DiscreteStateModel <: BranchModel end
abstract type SimulationModel <: BranchModel end #Simulation models typically can't propogate uncertainty, and aren't used for inference

abstract type StatePath end

#include("core/core.jl")
include("core/nodes/nodes.jl")
include("core/algorithms/algorithms.jl")
include("core/sim_tree.jl")
include("models/models.jl")
include("utils/utils.jl")

#Optional dependencies
function __init__()
    #Have FASTX installed if you want to read/write fasta files using our convenience functions
    @require FASTX = "c2308a5c-f048-11e8-3e8a-31650f418d12" include("utils/fasta_io.jl")

    #Have Compose and Colors installed if you want to draw trees
    @require Compose = "a81c6b42-2e10-5240-aca2-a61377ecd94b" begin
        @require Colors = "5ae59095-9a9b-59fe-a467-6f913c188581" begin
            using .Compose
            using .Colors
            include("viz/tree_compose.jl")
        end
    end

    #Have Phylo and Plots installed if you want to draw trees with Phylo
    @require Phylo = "aea672f4-3940-5932-aa44-993d1c3ff149" begin
        @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
            import .Phylo
            import .Plots
            include("viz/phylo_glue.jl")
        end
    end
end

export
    #Tree
    FelNode,
    BranchModel,
    Partition,

    #Partitions
    DiscretePartition,
    GaussianPartition,
    CustomDiscretePartition,
    NucleotidePartition,
    AminoAcidPartition,
    GappyNucleotidePartition,
    GappyAminoAcidPartition,
    CodonPartition,
    MultiSitePartition,

    #Models
    DiscreteStateModel,
    GaussianModel,
    BrownianMotion,
    DiagonalizedCTMC,
    GeneralCTMC,


    #core functions
    sim_tree,
    internal_message_init!,
    forward!,
    backward!,
    combine!,
    felsenstein!,
    felsenstein_down!,
    sample_down!,
    #endpoint_conditioned_sample_down!,
    log_likelihood!,
    marginal_state_dict,

    #tree functions
    ladderize,
    ladderize!,
    getleaflist,
    getnodelist,
    getnonleaflist,
    isleafnode,
    isrootnode,
    isinternalnode,
    isbranchnode,
    reroot!,
    nni_optim!,
    branchlength_optim!,

    #util functions
    one_hot_sample,
    scaled_prob_domain,
    golden_section_maximize,
    brents_method_minimize
    unit_transform,
    HKY85,
    P_from_diagonalized_Q,
    scale_cols_by_vec!,

    #things the user might overload
    copy_partition_to!,
    equilibrium_message,
    sample_partition!,
    obs2partition!,
    eq_freq,
    populate_tree!,
    populate_message!,
    partition2obs

end

#BM ToDo list:
#1. Flesh out the family of ancestral state algorithms, including:
#   a. Marginal ancestral state reconstruction
#   b. "Maximum joint probability" ancestral state reconstruction - requires "max" of partition.
#   - Is there a version of this that maintains multiple options whenever alternative states have non-negligible probability?
#   - In cases where a parent is ambiguous, but a child isn't, we don't want to duplicate the entire tree, so can we collapse back down?
#   - How would we represent this data structure. It seems like a DAG over trees, or something?
#   c. Sample ancestral states conditioned on the endpoint - requires "sample" of partition.
#These should all have consistent behavior in terms of what gets returned

#Personal reading list:
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5395464/
