
# Full API {#Full-API}

## Index {#Index}
- [`MolecularEvolution.AbstractUpdate`](#MolecularEvolution.AbstractUpdate)
- [`MolecularEvolution.AminoAcidPartition`](#MolecularEvolution.AminoAcidPartition)
- [`MolecularEvolution.BWMModel`](#MolecularEvolution.BWMModel-Union{Tuple{Vector{<:M}},%20Tuple{M}}%20where%20M<:DiscreteStateModel)
- [`MolecularEvolution.BWMModel`](#MolecularEvolution.BWMModel)
- [`MolecularEvolution.BranchlengthSampler`](#MolecularEvolution.BranchlengthSampler)
- [`MolecularEvolution.BrownianMotion`](#MolecularEvolution.BrownianMotion)
- [`MolecularEvolution.CATModel`](#MolecularEvolution.CATModel)
- [`MolecularEvolution.CATPartition`](#MolecularEvolution.CATPartition)
- [`MolecularEvolution.CladeStats`](#MolecularEvolution.CladeStats)
- [`MolecularEvolution.CodonPartition`](#MolecularEvolution.CodonPartition)
- [`MolecularEvolution.CovarionModel`](#MolecularEvolution.CovarionModel)
- [`MolecularEvolution.CovarionPartition`](#MolecularEvolution.CovarionPartition)
- [`MolecularEvolution.CustomDiscretePartition`](#MolecularEvolution.CustomDiscretePartition)
- [`MolecularEvolution.DiagonalizedCTMC`](#MolecularEvolution.DiagonalizedCTMC)
- [`MolecularEvolution.DiscretePartition`](#MolecularEvolution.DiscretePartition)
- [`MolecularEvolution.GappyAminoAcidPartition`](#MolecularEvolution.GappyAminoAcidPartition)
- [`MolecularEvolution.GappyNucleotidePartition`](#MolecularEvolution.GappyNucleotidePartition)
- [`MolecularEvolution.GaussianPartition`](#MolecularEvolution.GaussianPartition)
- [`MolecularEvolution.GeneralCTMC`](#MolecularEvolution.GeneralCTMC)
- [`MolecularEvolution.InterpolatedDiscreteModel`](#MolecularEvolution.InterpolatedDiscreteModel)
- [`MolecularEvolution.LazyDown`](#MolecularEvolution.LazyDown)
- [`MolecularEvolution.LazyPartition`](#MolecularEvolution.LazyPartition)
- [`MolecularEvolution.LazyUp`](#MolecularEvolution.LazyUp)
- [`MolecularEvolution.NucleotidePartition`](#MolecularEvolution.NucleotidePartition)
- [`MolecularEvolution.PModel`](#MolecularEvolution.PModel)
- [`MolecularEvolution.PiQ`](#MolecularEvolution.PiQ)
- [`MolecularEvolution.SWMModel`](#MolecularEvolution.SWMModel)
- [`MolecularEvolution.SWMPartition`](#MolecularEvolution.SWMPartition)
- [`MolecularEvolution.StandardUpdate`](#MolecularEvolution.StandardUpdate)
- [`MolecularEvolution.ZeroDriftOrnsteinUhlenbeck`](#MolecularEvolution.ZeroDriftOrnsteinUhlenbeck)
- [`Base.:==`](#Base.:==-Union{Tuple{T},%20Tuple{T,%20T}}%20where%20T<:AbstractTreeNode)
- [`MolecularEvolution.BayesUpdate`](#MolecularEvolution.BayesUpdate-Tuple{})
- [`MolecularEvolution.HB98_AAfit`](#MolecularEvolution.HB98_AAfit-Tuple{Any,%20Any,%20Any})
- [`MolecularEvolution.HB98_F61`](#MolecularEvolution.HB98_F61-Tuple{Any,%20Any,%20Any})
- [`MolecularEvolution.HIPSTR`](#MolecularEvolution.HIPSTR-Tuple{Vector{FelNode}})
- [`MolecularEvolution.HIPSTR`](#MolecularEvolution.HIPSTR)
- [`MolecularEvolution.MG94_F3x4`](#MolecularEvolution.MG94_F3x4-NTuple{4,%20Any})
- [`MolecularEvolution.MG94_F61`](#MolecularEvolution.MG94_F61-NTuple{4,%20Any})
- [`MolecularEvolution.MaxLikUpdate`](#MolecularEvolution.MaxLikUpdate-Tuple{})
- [`MolecularEvolution.SWM_prob_grid`](#MolecularEvolution.SWM_prob_grid-Union{Tuple{SWMPartition{PType}},%20Tuple{PType}}%20where%20PType<:MultiSitePartition)
- [`MolecularEvolution._mapreduce`](#MolecularEvolution._mapreduce-Union{Tuple{T},%20Tuple{AbstractTreeNode,%20T,%20Any,%20Any}}%20where%20T<:Function)
- [`MolecularEvolution.allocate!`](#MolecularEvolution.allocate!-Tuple{Any,%20Any})
- [`MolecularEvolution.backward!`](#MolecularEvolution.backward!-Tuple{DiscretePartition,%20DiscretePartition,%20MolecularEvolution.PMatrixModel,%20FelNode})
- [`MolecularEvolution.backward!`](#MolecularEvolution.backward!)
- [`MolecularEvolution.bfs_mapreduce`](#MolecularEvolution.bfs_mapreduce-Union{Tuple{T},%20Tuple{AbstractTreeNode,%20T,%20Any}}%20where%20T<:Function)
- [`MolecularEvolution.branchlength_optim!`](#MolecularEvolution.branchlength_optim!-Tuple)
- [`MolecularEvolution.branchlength_optim!`](#MolecularEvolution.branchlength_optim!)
- [`MolecularEvolution.branchlength_update!`](#MolecularEvolution.branchlength_update!-Tuple{MolecularEvolution.UnivariateModifier,%20FelNode,%20Any})
- [`MolecularEvolution.brents_method_minimize`](#MolecularEvolution.brents_method_minimize-Tuple{Any,%20Real,%20Real,%20Any,%20Real})
- [`MolecularEvolution.cascading_max_state_dict`](#MolecularEvolution.cascading_max_state_dict)
- [`MolecularEvolution.cascading_max_state_dict`](#MolecularEvolution.cascading_max_state_dict-Tuple{FelNode,%20Any})
- [`MolecularEvolution.char_proportions`](#MolecularEvolution.char_proportions-Tuple{Any,%20String})
- [`MolecularEvolution.collect_clades!`](#MolecularEvolution.collect_clades!-Tuple{FelNode,%20Dict{Tuple{UInt64,%20UInt64},%20MolecularEvolution.CladeStats},%20Dict{String,%20Int64}})
- [`MolecularEvolution.collect_leaf_dists`](#MolecularEvolution.collect_leaf_dists-Tuple{Vector{<:AbstractTreeNode}})
- [`MolecularEvolution.colored_seq_draw`](#MolecularEvolution.colored_seq_draw-Tuple{Any,%20Any,%20AbstractString})
- [`MolecularEvolution.combine!`](#MolecularEvolution.combine!)
- [`MolecularEvolution.combine!`](#MolecularEvolution.combine!-Tuple{DiscretePartition,%20DiscretePartition})
- [`MolecularEvolution.copy_tree`](#MolecularEvolution.copy_tree)
- [`MolecularEvolution.count_F3x4`](#MolecularEvolution.count_F3x4-Tuple{Array{String}})
- [`MolecularEvolution.count_F61`](#MolecularEvolution.count_F61-Tuple{Array{String}})
- [`MolecularEvolution.deepequals`](#MolecularEvolution.deepequals-Union{Tuple{T},%20Tuple{T,%20T}}%20where%20T<:AbstractTreeNode)
- [`MolecularEvolution.dfs_mapreduce`](#MolecularEvolution.dfs_mapreduce-Union{Tuple{T},%20Tuple{AbstractTreeNode,%20T,%20Any}}%20where%20T<:Function)
- [`MolecularEvolution.discrete_name_color_dict`](#MolecularEvolution.discrete_name_color_dict-Tuple{AbstractTreeNode,%20Any})
- [`MolecularEvolution.draw_example_tree`](#MolecularEvolution.draw_example_tree-Tuple{})
- [`MolecularEvolution.endpoint_conditioned_sample_state_dict`](#MolecularEvolution.endpoint_conditioned_sample_state_dict)
- [`MolecularEvolution.endpoint_conditioned_sample_state_dict`](#MolecularEvolution.endpoint_conditioned_sample_state_dict-Tuple{FelNode,%20Any})
- [`MolecularEvolution.expected_subs_per_site`](#MolecularEvolution.expected_subs_per_site-Tuple{Any,%20Any})
- [`MolecularEvolution.felsenstein!`](#MolecularEvolution.felsenstein!-Tuple{FelNode,%20Any})
- [`MolecularEvolution.felsenstein_down!`](#MolecularEvolution.felsenstein_down!-Tuple{FelNode,%20Any})
- [`MolecularEvolution.felsenstein_roundtrip!`](#MolecularEvolution.felsenstein_roundtrip!-Tuple{FelNode,%20Any})
- [`MolecularEvolution.forward!`](#MolecularEvolution.forward!-Tuple{DiscretePartition,%20DiscretePartition,%20MolecularEvolution.PMatrixModel,%20FelNode})
- [`MolecularEvolution.forward!`](#MolecularEvolution.forward!)
- [`MolecularEvolution.gappy_Q_from_symmetric_rate_matrix`](#MolecularEvolution.gappy_Q_from_symmetric_rate_matrix-Tuple{Any,%20Any,%20Any})
- [`MolecularEvolution.get_highlighter_legend`](#MolecularEvolution.get_highlighter_legend-Tuple{Any})
- [`MolecularEvolution.get_max_depth`](#MolecularEvolution.get_max_depth-Tuple{Any,%20Real})
- [`MolecularEvolution.get_phylo_tree`](#MolecularEvolution.get_phylo_tree)
- [`MolecularEvolution.get_phylo_tree`](#MolecularEvolution.get_phylo_tree-Tuple{FelNode})
- [`MolecularEvolution.golden_section_maximize`](#MolecularEvolution.golden_section_maximize-Tuple{Any,%20Real,%20Real,%20Any,%20Real})
- [`MolecularEvolution.hash_clade`](#MolecularEvolution.hash_clade-Tuple{BitSet})
- [`MolecularEvolution.highlight_seq_draw`](#MolecularEvolution.highlight_seq_draw-Tuple{Any,%20Any,%20AbstractString,%20Any,%20Any,%20Any})
- [`MolecularEvolution.highlighter_tree_draw`](#MolecularEvolution.highlighter_tree_draw-NTuple{4,%20Any})
- [`MolecularEvolution.internal_message_init!`](#MolecularEvolution.internal_message_init!-Tuple{FelNode,%20Partition})
- [`MolecularEvolution.internal_message_init!`](#MolecularEvolution.internal_message_init!-Tuple{FelNode,%20Vector{<:Partition}})
- [`MolecularEvolution.internal_nodes`](#MolecularEvolution.internal_nodes-Tuple{Any})
- [`MolecularEvolution.istreeconsistent`](#MolecularEvolution.istreeconsistent-Tuple{T}%20where%20T<:AbstractTreeNode)
- [`MolecularEvolution.ladder_tree_sim`](#MolecularEvolution.ladder_tree_sim-Tuple{Any})
- [`MolecularEvolution.lazyprep!`](#MolecularEvolution.lazyprep!-Tuple{FelNode,%20Vector{<:Partition}})
- [`MolecularEvolution.lazysort!`](#MolecularEvolution.lazysort!-Tuple{FelNode})
- [`MolecularEvolution.leaf_distmat`](#MolecularEvolution.leaf_distmat-Tuple{Any})
- [`MolecularEvolution.leaf_names`](#MolecularEvolution.leaf_names-Tuple{FelNode})
- [`MolecularEvolution.leaf_samples`](#MolecularEvolution.leaf_samples-Tuple{FelNode})
- [`MolecularEvolution.leaves`](#MolecularEvolution.leaves-Tuple{Any})
- [`MolecularEvolution.linear_scale`](#MolecularEvolution.linear_scale-NTuple{5,%20Any})
- [`MolecularEvolution.log_likelihood`](#MolecularEvolution.log_likelihood-Tuple{FelNode,%20BranchModel})
- [`MolecularEvolution.log_likelihood!`](#MolecularEvolution.log_likelihood!-Tuple{FelNode,%20Any})
- [`MolecularEvolution.longest_path`](#MolecularEvolution.longest_path-Tuple{FelNode})
- [`MolecularEvolution.marginal_state_dict`](#MolecularEvolution.marginal_state_dict)
- [`MolecularEvolution.marginal_state_dict`](#MolecularEvolution.marginal_state_dict-Tuple{FelNode,%20Any})
- [`MolecularEvolution.matrix_for_display`](#MolecularEvolution.matrix_for_display-Tuple{Any,%20Any})
- [`MolecularEvolution.metropolis_sample`](#MolecularEvolution.metropolis_sample-Tuple{FelNode,%20Vector{<:BranchModel},%20Any})
- [`MolecularEvolution.metropolis_sample`](#MolecularEvolution.metropolis_sample-Tuple{AbstractUpdate,%20FelNode,%20Any,%20Any})
- [`MolecularEvolution.metropolis_step`](#MolecularEvolution.metropolis_step-Tuple{Function,%20Any,%20Any})
- [`MolecularEvolution.midpoint`](#MolecularEvolution.midpoint-Tuple{FelNode})
- [`MolecularEvolution.mix`](#MolecularEvolution.mix-Union{Tuple{SWMPartition{PType}},%20Tuple{PType}}%20where%20PType<:DiscretePartition)
- [`MolecularEvolution.name2node_dict`](#MolecularEvolution.name2node_dict-Tuple{AbstractTreeNode})
- [`MolecularEvolution.newick`](#MolecularEvolution.newick)
- [`MolecularEvolution.newick`](#MolecularEvolution.newick-Tuple{AbstractTreeNode})
- [`MolecularEvolution.nni_optim!`](#MolecularEvolution.nni_optim!)
- [`MolecularEvolution.nni_optim!`](#MolecularEvolution.nni_optim!-Tuple)
- [`MolecularEvolution.nni_update!`](#MolecularEvolution.nni_update!-Tuple{Function,%20FelNode,%20Any})
- [`MolecularEvolution.node_distances`](#MolecularEvolution.node_distances-Tuple{AbstractTreeNode})
- [`MolecularEvolution.node_names`](#MolecularEvolution.node_names-Tuple{FelNode})
- [`MolecularEvolution.node_samples`](#MolecularEvolution.node_samples-Tuple{FelNode})
- [`MolecularEvolution.nodes`](#MolecularEvolution.nodes-Tuple{Any})
- [`MolecularEvolution.nonreversibleQ`](#MolecularEvolution.nonreversibleQ-Tuple{Any})
- [`MolecularEvolution.parent_list`](#MolecularEvolution.parent_list-Tuple{FelNode})
- [`MolecularEvolution.partition2obs`](#MolecularEvolution.partition2obs-Tuple{DiscretePartition,%20String})
- [`MolecularEvolution.partition2obs`](#MolecularEvolution.partition2obs)
- [`MolecularEvolution.plot_multiple_trees`](#MolecularEvolution.plot_multiple_trees)
- [`MolecularEvolution.plot_multiple_trees`](#MolecularEvolution.plot_multiple_trees-Tuple{Any,%20Any})
- [`MolecularEvolution.populate_tree!`](#MolecularEvolution.populate_tree!-Tuple{FelNode,%20Partition,%20Any,%20Any})
- [`MolecularEvolution.populate_tree!`](#MolecularEvolution.populate_tree!)
- [`MolecularEvolution.promote_internal`](#MolecularEvolution.promote_internal-Tuple{FelNode})
- [`MolecularEvolution.quadratic_CI`](#MolecularEvolution.quadratic_CI-Tuple{Vector,%20Vector})
- [`MolecularEvolution.quadratic_CI`](#MolecularEvolution.quadratic_CI-Tuple{Function,%20Vector,%20Int64})
- [`MolecularEvolution.read_fasta`](#MolecularEvolution.read_fasta)
- [`MolecularEvolution.read_fasta`](#MolecularEvolution.read_fasta-Tuple{String})
- [`MolecularEvolution.read_newick_tree`](#MolecularEvolution.read_newick_tree)
- [`MolecularEvolution.read_newick_tree`](#MolecularEvolution.read_newick_tree-Tuple{String})
- [`MolecularEvolution.refresh!`](#MolecularEvolution.refresh!-Tuple{FelNode,%20Any})
- [`MolecularEvolution.reversibleQ`](#MolecularEvolution.reversibleQ)
- [`MolecularEvolution.reversibleQ`](#MolecularEvolution.reversibleQ-Tuple{Any,%20Any})
- [`MolecularEvolution.root2tip_distances`](#MolecularEvolution.root2tip_distances-Tuple{AbstractTreeNode})
- [`MolecularEvolution.root_optim!`](#MolecularEvolution.root_optim!)
- [`MolecularEvolution.root_optim!`](#MolecularEvolution.root_optim!-Tuple{FelNode,%20Any})
- [`MolecularEvolution.root_update!`](#MolecularEvolution.root_update!-Tuple{RootUpdate,%20FelNode,%20Any})
- [`MolecularEvolution.sample_down!`](#MolecularEvolution.sample_down!)
- [`MolecularEvolution.sample_down!`](#MolecularEvolution.sample_down!-Tuple{FelNode,%20Any,%20Any})
- [`MolecularEvolution.sample_from_message!`](#MolecularEvolution.sample_from_message!-Tuple{Vector{<:Partition}})
- [`MolecularEvolution.savefig_tweakSVG`](#MolecularEvolution.savefig_tweakSVG-Tuple{Any,%20Context})
- [`MolecularEvolution.savefig_tweakSVG`](#MolecularEvolution.savefig_tweakSVG)
- [`MolecularEvolution.savefig_tweakSVG`](#MolecularEvolution.savefig_tweakSVG-Tuple{Any,%20Plots.Plot})
- [`MolecularEvolution.shortest_path_between_nodes`](#MolecularEvolution.shortest_path_between_nodes-Tuple{FelNode,%20FelNode})
- [`MolecularEvolution.sibling_inds`](#MolecularEvolution.sibling_inds-Tuple{AbstractTreeNode})
- [`MolecularEvolution.siblings`](#MolecularEvolution.siblings-Tuple{AbstractTreeNode})
- [`MolecularEvolution.sim_tree`](#MolecularEvolution.sim_tree-Tuple{})
- [`MolecularEvolution.sim_tree`](#MolecularEvolution.sim_tree-Tuple{Int64,%20Any,%20Any})
- [`MolecularEvolution.sim_tree`](#MolecularEvolution.sim_tree)
- [`MolecularEvolution.simple_radial_tree_plot`](#MolecularEvolution.simple_radial_tree_plot-Tuple{FelNode})
- [`MolecularEvolution.simple_tree_draw`](#MolecularEvolution.simple_tree_draw-Tuple{FelNode})
- [`MolecularEvolution.standard_tree_sim`](#MolecularEvolution.standard_tree_sim-Tuple{Any})
- [`MolecularEvolution.total_LL`](#MolecularEvolution.total_LL-Tuple{Partition})
- [`MolecularEvolution.tree2distances`](#MolecularEvolution.tree2distances-Tuple{AbstractTreeNode})
- [`MolecularEvolution.tree2shared_branch_lengths`](#MolecularEvolution.tree2shared_branch_lengths-Tuple{AbstractTreeNode})
- [`MolecularEvolution.tree_draw`](#MolecularEvolution.tree_draw)
- [`MolecularEvolution.tree_draw`](#MolecularEvolution.tree_draw-Tuple{FelNode})
- [`MolecularEvolution.tree_polish!`](#MolecularEvolution.tree_polish!)
- [`MolecularEvolution.tree_polish!`](#MolecularEvolution.tree_polish!-Tuple{Any,%20Any})
- [`MolecularEvolution.unc2probvec`](#MolecularEvolution.unc2probvec)
- [`MolecularEvolution.unc2probvec`](#MolecularEvolution.unc2probvec-Tuple{Any})
- [`MolecularEvolution.univariate_maximize`](#MolecularEvolution.univariate_maximize-Tuple{Any,%20Real,%20Real,%20Any,%20GoldenSectionOpt,%20Real})
- [`MolecularEvolution.univariate_maximize`](#MolecularEvolution.univariate_maximize-Tuple{Any,%20Real,%20Real,%20Any,%20BrentsMethodOpt,%20Real})
- [`MolecularEvolution.univariate_sampler`](#MolecularEvolution.univariate_sampler-Tuple{Any,%20BranchlengthSampler,%20Any})
- [`MolecularEvolution.values_from_phylo_tree`](#MolecularEvolution.values_from_phylo_tree-Tuple{Any,%20Any})
- [`MolecularEvolution.values_from_phylo_tree`](#MolecularEvolution.values_from_phylo_tree)
- [`MolecularEvolution.weightEM`](#MolecularEvolution.weightEM-Tuple{Matrix{Float64},%20Any})
- [`MolecularEvolution.write_fasta`](#MolecularEvolution.write_fasta)
- [`MolecularEvolution.write_fasta`](#MolecularEvolution.write_fasta-Tuple{String,%20Vector{String}})
- [`MolecularEvolution.write_nexus`](#MolecularEvolution.write_nexus)
- [`MolecularEvolution.write_nexus`](#MolecularEvolution.write_nexus-Tuple{String,%20FelNode})


## Docstrings {#Docstrings}
<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.AbstractUpdate' href='#MolecularEvolution.AbstractUpdate'><span class="jlbinding">MolecularEvolution.AbstractUpdate</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



**Summary**

`abstract type AbstractUpdate <: Function`

A callable type that typically takes `(tree::FelNode, models; partition_list=1:length(tree.message))`, updates `tree` and `models`, and returns the updated `tree` and `models`.

**Example**

Define a new subtype, where `foo` and `bar` are arbitrary updating functions

```julia
struct MyUpdate <: AbstractUpdate end

function (update::MyUpdate)(tree::FelNode, models; partition_list=1:length(tree.message))
    tree, models = foo(tree, models, partition_list=partition_list)
    tree, models = BayesUpdate(nni=0)(tree, models, partition_list=partition_list)
    tree, models = bar(tree, models, partition_list=partition_list)
    return tree, models
end
```


See also: [`StandardUpdate`](/api#MolecularEvolution.StandardUpdate)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/algorithms/AbstractUpdate.jl#L1-L20" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.AminoAcidPartition' href='#MolecularEvolution.AminoAcidPartition'><span class="jlbinding">MolecularEvolution.AminoAcidPartition</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
mutable struct AminoAcidPartition <: DiscretePartition
```


**Constructors**

```
AminoAcidPartition(sites)
AminoAcidPartition(freq_vec::Vector{Float64}, sites::Int64)
AminoAcidPartition(state, states, sites, scaling)
```


**Description**

A `DiscretePartition` for amino acid sequences with 20 states representing the standard amino acids.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/discrete_models/discrete_partitions.jl#L125-L134" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.BWMModel' href='#MolecularEvolution.BWMModel'><span class="jlbinding">MolecularEvolution.BWMModel</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
mutable struct BWMModel{M} <: DiscreteStateModel where M <: DiscreteStateModel
```


**Fields**
- `models::Vector{<:M}`: A vector of models.
  
- `weights::Vector{Float64}`: A vector of weights.
  

**Description**

Branch-wise mixture model.

::: tip Note

`forward!` and `backward!` are currently only defined for `M<:PMatrixModel`.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/compound_models/bwm.jl#L2-L13" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.BWMModel-Union{Tuple{Vector{<:M}}, Tuple{M}} where M<:DiscreteStateModel' href='#MolecularEvolution.BWMModel-Union{Tuple{Vector{<:M}}, Tuple{M}} where M<:DiscreteStateModel'><span class="jlbinding">MolecularEvolution.BWMModel</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
BWMModel{M}(models::Vector{<:M}) where M <: DiscreteStateModel
```


Convenience constructor where the weights are uniform.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/compound_models/bwm.jl#L19-L23" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.BranchlengthSampler' href='#MolecularEvolution.BranchlengthSampler'><span class="jlbinding">MolecularEvolution.BranchlengthSampler</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
BranchlengthSampler

A type that allows you to specify a additive proposal function in the log domain and a prior distrubution over the log of the branchlengths. It also stores `acc_ratio` which is a tuple of `(ratio, total, #acceptances)`, where `ratio::Float64` is the acceptance ratio, `total::Int64` is the total number of proposals, and `#acceptances::Int64` is the number of acceptances.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/utils/simple_sample.jl#L6-L10" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.BrownianMotion' href='#MolecularEvolution.BrownianMotion'><span class="jlbinding">MolecularEvolution.BrownianMotion</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
BrownianMotion(mean_drift::Float64, var_drift::Float64)
```


A 1D continuous Brownian motion model with mean drift and variance drift.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/continuous_models/brownian_motion.jl#L2-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.CATModel' href='#MolecularEvolution.CATModel'><span class="jlbinding">MolecularEvolution.CATModel</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
CATModel(models::Vector{<:BranchModel})
```


CAT is something where you split the sites up, and assign each site to a different model (whose &quot;data&quot; gets stored in a contiguous block of memory).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/compound_models/cat.jl#L2-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.CATPartition' href='#MolecularEvolution.CATPartition'><span class="jlbinding">MolecularEvolution.CATPartition</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



**Constructors**

```julia
CATPartition(part_inds::Vector{Vector{Int}})
CATPartition(part_inds::Vector{Vector{Int}}, parts::Vector{PType})
```


**Description**

A partition for the [`CATModel`](/api#MolecularEvolution.CATModel).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/compound_models/cat.jl#L15-L23" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.CladeStats' href='#MolecularEvolution.CladeStats'><span class="jlbinding">MolecularEvolution.CladeStats</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Store statistics about a clade: its frequency and observed child pairs.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/utils/HIPSTR.jl#L164-L166" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.CodonPartition' href='#MolecularEvolution.CodonPartition'><span class="jlbinding">MolecularEvolution.CodonPartition</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
mutable struct CodonPartition <: DiscretePartition
```


**Constructors**

```
CodonPartition(sites; code = universal_code)
CodonPartition(state, states, sites, scaling; code = universal_code)
```


**Description**

A `DiscretePartition` where every state represents a codon. Can be customized to use different genetic codes.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/discrete_models/codon_models.jl#L349-L357" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.CovarionModel' href='#MolecularEvolution.CovarionModel'><span class="jlbinding">MolecularEvolution.CovarionModel</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



**Constructors**

```julia
CovarionModel(models::Vector{<:DiscreteStateModel}, inter_transition_rates::Matrix{Float64})
CovarionModel(models::Vector{<:DiscreteStateModel}, inter_transition_rate::Float64)
```


**Description**

The covarion model.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/compound_models/covarion.jl#L2-L10" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.CovarionPartition' href='#MolecularEvolution.CovarionPartition'><span class="jlbinding">MolecularEvolution.CovarionPartition</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
CovarionPartition(states,sites,models,t)
```


A partition for the [`CovarionModel`](/api#MolecularEvolution.CovarionModel).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/compound_models/covarion.jl#L47-L51" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.CustomDiscretePartition' href='#MolecularEvolution.CustomDiscretePartition'><span class="jlbinding">MolecularEvolution.CustomDiscretePartition</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
mutable struct CustomDiscretePartition <: DiscretePartition
```


**Constructors**

```
CustomDiscretePartition(states, sites)
CustomDiscretePartition(freq_vec::Vector{Float64}, sites::Int64)
CustomDiscretePartition(state, states, sites, scaling)
```


**Description**

A general-purpose `DiscretePartition` that can handle any number of states. Useful for custom discrete state spaces that don&#39;t fit into the predefined partition types.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/discrete_models/discrete_partitions.jl#L35-L45" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.DiagonalizedCTMC' href='#MolecularEvolution.DiagonalizedCTMC'><span class="jlbinding">MolecularEvolution.DiagonalizedCTMC</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



**Constructors**

```julia
DiagonalizedCTMC(Q::Array{Float64,2})
DiagonalizedCTMC(Qpre::Array{Float64,2}, pi::Vector{Float64})
```


**Description**

Takes in a Q matrix (which can be multiplied onto row-wise by `pi`) and diagonalizes it. When computing $e^{Q t}$ (for different $t$s), we now only need to exponentiate the eigenvalues, and perform the two change-of-basis matrix multiplications.

::: warning Warning

Construction fails if `Q` has complex eigenvalues.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/discrete_models/pmatrix_models/DiagonalizedCTMC.jl#L5-L17" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.DiscretePartition' href='#MolecularEvolution.DiscretePartition'><span class="jlbinding">MolecularEvolution.DiscretePartition</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
abstract type DiscretePartition <: MultiSitePartition end
```


Represents a `states`-by-`sites` matrix of probabilities. The following fields are loosely required:
- `state`: A matrix of probabilities that are site-wise normalized.
  
- `states`: The number of states.
  
- `sites`: The number of sites.
  
- `scaling`: A vector of log-domain probability scaling factors, one per site.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/discrete_models/discrete_partitions.jl#L5-L13" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.GappyAminoAcidPartition' href='#MolecularEvolution.GappyAminoAcidPartition'><span class="jlbinding">MolecularEvolution.GappyAminoAcidPartition</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
mutable struct GappyAminoAcidPartition <: DiscretePartition
```


**Constructors**

```
GappyAminoAcidPartition(sites)
GappyAminoAcidPartition(freq_vec::Vector{Float64}, sites::Int64)
GappyAminoAcidPartition(state, states, sites, scaling)
```


**Description**

A `DiscretePartition` for amino acid sequences with 21 states (20 standard amino acids plus &#39;-&#39; for gaps).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/discrete_models/discrete_partitions.jl#L155-L164" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.GappyNucleotidePartition' href='#MolecularEvolution.GappyNucleotidePartition'><span class="jlbinding">MolecularEvolution.GappyNucleotidePartition</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
mutable struct GappyNucleotidePartition <: DiscretePartition
```


**Constructors**

```
GappyNucleotidePartition(sites)
GappyNucleotidePartition(freq_vec::Vector{Float64}, sites::Int64)
GappyNucleotidePartition(state, states, sites, scaling)
```


**Description**

A `DiscretePartition` for nucleotide sequences with 5 states (A, C, G, T, -), where &#39;-&#39; represents a gap.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/discrete_models/discrete_partitions.jl#L95-L104" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.GaussianPartition' href='#MolecularEvolution.GaussianPartition'><span class="jlbinding">MolecularEvolution.GaussianPartition</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



**Constructors**

```julia
GaussianPartition(mean::Float64, var::Float64, norm_const::Float64)
GaussianPartition(mean::Float64, var::Float64) # norm_const defaults to 0.0
GaussianPartition() # mean, var, norm_const default to 0.0, 1.0, 0.0 respectively
```


**Description**

A partition representing a (not necessarily normalized) Gaussian distribution. `norm_const` is the log-domain normalization constant.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/continuous_models/gaussian_partition.jl#L3-L12" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.GeneralCTMC' href='#MolecularEvolution.GeneralCTMC'><span class="jlbinding">MolecularEvolution.GeneralCTMC</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
function GeneralCTMC(Q::Array{Float64,2})
```


Wraps a Q matrix and will compute $e^{Q t}$ (for different $t$s) naively by `exp`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/discrete_models/pmatrix_models/GeneralCTMC.jl#L1-L5" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.InterpolatedDiscreteModel' href='#MolecularEvolution.InterpolatedDiscreteModel'><span class="jlbinding">MolecularEvolution.InterpolatedDiscreteModel</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



**Constructors**

```julia
InterpolatedDiscreteModel(siz::Int64, generator, tvec::Vector{Float64})
InterpolatedDiscreteModel(Pvec::Array{Float64,3}, tvec::Vector{Float64})
```


`generator` is a function that takes a time value `t` and returns a P matrix.

**Description**

Stores a number (`siz`) of P matrices, and the time values to which they correspond. For a requested t, the returned P matrix is (element-wise linearly) interpolated between it&#39;s two neighbours.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/discrete_models/interpolated_discrete_model.jl#L5-L15" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.LazyDown' href='#MolecularEvolution.LazyDown'><span class="jlbinding">MolecularEvolution.LazyDown</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



**Constructors**

```
LazyDown(stores_obs)
LazyDown() = LazyDown(x::FelNode -> true)
```


**Description**

Indicate that we want to do a downward pass, e.g. `sample_down!`. The function passed to the constructor takes a `node::FelNode` as input and returns a `Bool` that decides if `node` stores its observations.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/lazy_models/lazy_partition.jl#L51-L59" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.LazyPartition' href='#MolecularEvolution.LazyPartition'><span class="jlbinding">MolecularEvolution.LazyPartition</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



**Constructor**

```
LazyPartition{PType}()
```


Initialize an empty `LazyPartition` that is meant for wrapping a partition of type `PType`.

**Description**

With this data structure, you can wrap a partition of choice.  The idea is that in some message passing algorithms, there is only a wave of partitions which need to actualize.  For instance, a wave following a root-leaf path, or a depth-first traversal. In which case, we can be more economical with our memory consumption. With a worst case memory complexity of O(log(n)), where n is the number of nodes, functionality is provided for:
- `log_likelihood!`
  
- `felsenstein!`
  
- `sample_down!`
  

::: tip Note

For successive `felsenstein!` calls, we need to extract the information at the root somehow after each call. This can be done with e.g. `total_LL` or `site_LLs`.

:::

**Further requirements**

Suppose you want to wrap a partition of `PType` with `LazyPartition`:
- If you&#39;re calling `log_likelihood!` and `felsenstein!`:
  - `obs2partition!(partition::PType, obs)` that transforms an observation to a partition.
    
  
- If you&#39;re calling `sample_down!`:
  - `partition2obs(partition::PType)` that returns the most likely state from a partition, inverts `obs2partition!`.
    
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/lazy_models/lazy_partition.jl#L2-L25" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.LazyUp' href='#MolecularEvolution.LazyUp'><span class="jlbinding">MolecularEvolution.LazyUp</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



**Constructor**

```
LazyUp()
```


**Description**

Indicate that we want to do an upward pass, e.g. `felsenstein!`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/lazy_models/lazy_partition.jl#L42-L48" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.NucleotidePartition' href='#MolecularEvolution.NucleotidePartition'><span class="jlbinding">MolecularEvolution.NucleotidePartition</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
mutable struct NucleotidePartition <: DiscretePartition
```


**Constructors**

```
NucleotidePartition(sites)
NucleotidePartition(freq_vec::Vector{Float64}, sites::Int64)
NucleotidePartition(state, states, sites, scaling)
```


**Description**

A `DiscretePartition` specifically for nucleotide sequences with 4 states (A, C, G, T).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/discrete_models/discrete_partitions.jl#L65-L74" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.PModel' href='#MolecularEvolution.PModel'><span class="jlbinding">MolecularEvolution.PModel</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
PModel(P::Array{Float64,2})
```


A P matrix.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/discrete_models/pmatrix_models/PModel.jl#L2-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.PiQ' href='#MolecularEvolution.PiQ'><span class="jlbinding">MolecularEvolution.PiQ</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



**Constructors**

```julia
PiQ(r::Float64,pi::Vector{Float64}; normalize=false)
PiQ(pi::Vector{Float64}; normalize=false)
```


**Description**

The F81 substitution model, but for general dimensions. https://www.diva-portal.org/smash/get/diva2:1878793/FULLTEXT01.pdf


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/discrete_models/PiQ.jl#L4-L13" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.SWMModel' href='#MolecularEvolution.SWMModel'><span class="jlbinding">MolecularEvolution.SWMModel</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



**Constructors**

```julia
SWMModel(models::Vector{<:BranchModel})
SWMModel(model::M, rs::Vector{Float64}) where {M <: BranchModel}
```


**Description**

A site-wise mixture model, for site-to-site &quot;random effects&quot; rate variation.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/compound_models/swm.jl#L17-L25" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.SWMPartition' href='#MolecularEvolution.SWMPartition'><span class="jlbinding">MolecularEvolution.SWMPartition</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



**Constructors**

```julia
SWMPartition(parts::Vector{PType}) where {PType <: MultiSitePartition}
SWMPartition(part::PType, n_parts::Int) where {PType <: MultiSitePartition}
SWMPartition(parts::Vector{PType}, weights::Vector{Float64}, sites::Int, states::Int, models::Int) where {PType <: MultiSitePartition}
```


**Description**

A site-wise mixture partition for the [`SWMModel`](/api#MolecularEvolution.SWMModel).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/compound_models/swm.jl#L46-L55" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.StandardUpdate' href='#MolecularEvolution.StandardUpdate'><span class="jlbinding">MolecularEvolution.StandardUpdate</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



**Summary**

`struct StandardUpdate <: AbstractUpdate`

A standard update can be a family of calls to [`nni_update!`](/api#MolecularEvolution.nni_update!-Tuple{Function,%20FelNode,%20Any}), [`branchlength_update!`](/api#MolecularEvolution.branchlength_update!-Tuple{MolecularEvolution.UnivariateModifier,%20FelNode,%20Any}), [`root_update!`](/api#MolecularEvolution.root_update!-Tuple{RootUpdate,%20FelNode,%20Any}), and model updates.

**Constructor**

```
StandardUpdate(
    nni::Int,
    branchlength::Int,
    root::Int,
    models::Int,
    refresh::Bool,
    nni_selection::Function,
    branchlength_modifier::UnivariateModifier,
    root_update::RootUpdate,
    models_update::ModelsUpdate
)
```


**Arguments**
- `nni::Int`: the number of times to update the tree by `nni_update!`
  
- `branchlength::Int`: the number of times to update the tree by `branchlength_update!`
  
- `root::Int`: the number of times to update the tree by `root_update!`
  
- `models::Int`: the number of times to update the model
  
- `refresh::Bool`: whether to refresh the messages in tree between update operations to ensure message consistency
  
- `nni_selection::Function`: the function that selects between nni configurations
  
- `branchlength_modifier::UnivariateModifier`: the modifier to update a branchlength by `branchlength_update!`
  
- `root_update::RootUpdate`: updates the root by `root_update!`
  
- `models_update::ModelsUpdate`: updates the model parameters
  

See also: [`BayesUpdate`](/api#MolecularEvolution.BayesUpdate-Tuple{}), [`MaxLikUpdate`](/api#MolecularEvolution.MaxLikUpdate-Tuple{})


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/algorithms/AbstractUpdate.jl#L23-L53" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.ZeroDriftOrnsteinUhlenbeck' href='#MolecularEvolution.ZeroDriftOrnsteinUhlenbeck'><span class="jlbinding">MolecularEvolution.ZeroDriftOrnsteinUhlenbeck</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
function ZeroDriftOrnsteinUhlenbeck(
    mean::Float64,
    volatility::Float64,
    reversion::Float64,
)
```


A 1D continuous Ornstein-Uhlenbeck process with mean, volatility, and reversion.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/continuous_models/ornstein_uhlenbeck.jl#L3-L11" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Base.:==-Union{Tuple{T}, Tuple{T, T}} where T<:AbstractTreeNode' href='#Base.:==-Union{Tuple{T}, Tuple{T, T}} where T<:AbstractTreeNode'><span class="jlbinding">Base.:==</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
==(t1, t2)
Defaults to pointer equality
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/nodes/AbstractTreeNode.jl#L78-L81" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.BayesUpdate-Tuple{}' href='#MolecularEvolution.BayesUpdate-Tuple{}'><span class="jlbinding">MolecularEvolution.BayesUpdate</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
BayesUpdate(;
    nni = 1,
    branchlength = 1,
    root = 0,
    models = 0,
    refresh = false,
    branchlength_sampler::UnivariateSampler = BranchlengthSampler(
        Normal(0, 2),
        Normal(-1, 1),
    ),
    root_sampler::RootSample = StandardRootSample(1.0, 1),
    models_sampler::ModelsUpdate = StandardModelsUpdate()
)
```


Convenience constructor for [`StandardUpdate`](/api#MolecularEvolution.StandardUpdate). The `nni_selection` is fixed to `softmax_sampler`.  This constructor provides Bayesian updates by sampling from the posterior distribution.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/algorithms/AbstractUpdate.jl#L66-L83" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.HB98_AAfit-Tuple{Any, Any, Any}' href='#MolecularEvolution.HB98_AAfit-Tuple{Any, Any, Any}'><span class="jlbinding">MolecularEvolution.HB98_AAfit</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
function HB98_AAfit(alpha, nuc_matrix, AA_fitness; genetic_code = universal_code)
```


Compute the (Halpern and Bruno 1998) codon substitution Q matrix parameterized by `alpha`, a nucleotide substitution matrix `nuc_matrix`, and direct amino acid fitness `AA_fitness` vector.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/discrete_models/codon_models.jl#L304-L309" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.HB98_F61-Tuple{Any, Any, Any}' href='#MolecularEvolution.HB98_F61-Tuple{Any, Any, Any}'><span class="jlbinding">MolecularEvolution.HB98_F61</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
function HB98_F61(alpha, nuc_matrix, F61; genetic_code = universal_code)
```


Compute the (Halpern and Bruno 1998) HB98-F61 codon substitution Q matrix parameterized by `alpha`, a nucleotide substitution matrix `nuc_matrix`, and the `F61` codon frequency estimator (see `?MolecularEvolution.count_F61`).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/discrete_models/codon_models.jl#L264-L269" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.HIPSTR-Tuple{Vector{FelNode}}' href='#MolecularEvolution.HIPSTR-Tuple{Vector{FelNode}}'><span class="jlbinding">MolecularEvolution.HIPSTR</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
HIPSTR(trees::Vector{FelNode}; set_branchlengths = true)
```


Construct a Highest Independent Posterior Subtree Reconstruction (HIPSTR) tree from a collection of trees.

Returns a single FelNode representing the HIPSTR consensus tree.

If `set_branchlengths = true`, the branch length of a node in the HIPSTR tree will be set to the mean branch length of all nodes from the input trees that have the same clade. (By the same clade, we mean that the set of leaves below the node is the same.) Otherwise, the root branch length is 0.0 and the rest 1.0.

Source: https://www.biorxiv.org/content/10.1101/2024.12.08.627395v1.full.pdf


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/utils/HIPSTR.jl#L2-L13" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.MG94_F3x4-NTuple{4, Any}' href='#MolecularEvolution.MG94_F3x4-NTuple{4, Any}'><span class="jlbinding">MolecularEvolution.MG94_F3x4</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
function MG94_F3x4(alpha, beta, nuc_matrix, F3x4; genetic_code = universal_code)
```


Compute the (Muse and Gaut 1994) MG94-F3x4 codon substitution Q matrix parameterized by `alpha`, `beta`, a nucleotide substitution matrix `nuc_matrix`, and `F3x4` (see `?MolecularEvolution.count_F3x4`).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/discrete_models/codon_models.jl#L202-L207" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.MG94_F61-NTuple{4, Any}' href='#MolecularEvolution.MG94_F61-NTuple{4, Any}'><span class="jlbinding">MolecularEvolution.MG94_F61</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
function MG94_F61(alpha, beta, nuc_matrix, F61; genetic_code = universal_code)
```


Compute the (Muse and Gaut 1994) MG94-F61 codon substitution Q matrix parameterized by `alpha`, `beta`, a nucleotide substitution matrix `nuc_matrix`, and the `F61` codon frequency estimator (see `?MolecularEvolution.count_F61`).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/discrete_models/codon_models.jl#L243-L248" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.MaxLikUpdate-Tuple{}' href='#MolecularEvolution.MaxLikUpdate-Tuple{}'><span class="jlbinding">MolecularEvolution.MaxLikUpdate</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
MaxLikUpdate(;
    nni = 1,
    branchlength = 1,
    root = 0,
    models = 0,
    refresh = false,
    branchlength_optimizer::UnivariateOpt = GoldenSectionOpt(),
    root_optimizer = StandardRootOpt(10),
    models_optimizer::ModelUpdate = StandardModelUpdate()
)
```


Convenience constructor for [`StandardUpdate`](/api#MolecularEvolution.StandardUpdate). The `nni_selection` is fixed to `argmax`. This constructor provides Maximum Likelihood updates by optimizing parameters.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/algorithms/AbstractUpdate.jl#L108-L122" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.SWM_prob_grid-Union{Tuple{SWMPartition{PType}}, Tuple{PType}} where PType<:MultiSitePartition' href='#MolecularEvolution.SWM_prob_grid-Union{Tuple{SWMPartition{PType}}, Tuple{PType}} where PType<:MultiSitePartition'><span class="jlbinding">MolecularEvolution.SWM_prob_grid</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
SWM_prob_grid(part::SWMPartition{PType}) where {PType <: MultiSitePartition}
```


Returns a matrix of probabilities for each site, for each model (in the probability domain - not logged!) as well as the log probability offsets


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/compound_models/swm.jl#L150-L154" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution._mapreduce-Union{Tuple{T}, Tuple{AbstractTreeNode, T, Any, Any}} where T<:Function' href='#MolecularEvolution._mapreduce-Union{Tuple{T}, Tuple{AbstractTreeNode, T, Any, Any}} where T<:Function'><span class="jlbinding">MolecularEvolution._mapreduce</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Internal function. Helper for bfs_mapreduce and dfs_mapreduce


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/nodes/base_tree_utils.jl#L1" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.allocate!-Tuple{Any, Any}' href='#MolecularEvolution.allocate!-Tuple{Any, Any}'><span class="jlbinding">MolecularEvolution.allocate!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
allocate!(tree, partition_or_message)
```


Allocates initial messages for all nodes in the tree, copying the passed-in message template. If passed a partition, then this will assume the message template is a vector containing just that partition.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/simple_interface.jl#L38-L43" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.backward!-Tuple{DiscretePartition, DiscretePartition, MolecularEvolution.PMatrixModel, FelNode}' href='#MolecularEvolution.backward!-Tuple{DiscretePartition, DiscretePartition, MolecularEvolution.PMatrixModel, FelNode}'><span class="jlbinding">MolecularEvolution.backward!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
backward!(dest::Partition, source::Partition, model::BranchModel, node::FelNode)
```


Propagate the source partition backwards along the branch to the destination partition, under the model. Note: You should overload this for your own BranchModel types.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/discrete_models/pmatrix_models/pmatrix_models.jl#L5-L10" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.bfs_mapreduce-Union{Tuple{T}, Tuple{AbstractTreeNode, T, Any}} where T<:Function' href='#MolecularEvolution.bfs_mapreduce-Union{Tuple{T}, Tuple{AbstractTreeNode, T, Any}} where T<:Function'><span class="jlbinding">MolecularEvolution.bfs_mapreduce</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Performs a BFS map-reduce over the tree, starting at a given node For each node, map_reduce is called as:    map_reduce(curr_node::FelNode, prev_node::FelNode, aggregator) where prev_node is the previous node visited on the path from the start node to the current node It is expected to update the aggregator, and not return anything.

Not exactly conventional map-reduce, as map-reduce calls may rely on state in the aggregator added by map-reduce calls on other nodes visited earlier.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/nodes/base_tree_utils.jl#L27-L36" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.branchlength_optim!-Tuple' href='#MolecularEvolution.branchlength_optim!-Tuple'><span class="jlbinding">MolecularEvolution.branchlength_optim!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
branchlength_optim!(tree::FelNode, models;  <keyword arguments>)
```


Uses golden section search, or optionally Brent&#39;s method, to optimize all branches recursively, maintaining the integrity of the messages. Requires felsenstein!() to have been run first. models can either be a single model (if the messages on the tree contain just one Partition) or an array of models, if the messages have &gt;1 Partition, or  a function that takes a node, and returns a Vector{&lt;:BranchModel} if you need the models to vary from one branch to another.

**Keyword Arguments**
- `partition_list=nothing`: (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over (but you probably want to optimize branch lengths with all models, the default option).
  
- `tol=1e-5`: absolute tolerance for the `bl_optimizer`.
  
- `bl_optimizer::UnivariateModifier=GoldenSectionOpt()`: the algorithm used to optimize the log likelihood of a branch length. In addition to golden section search, Brent&#39;s method can be used by setting `bl_optimizer=BrentsMethodOpt()`.
  
- `sort_tree=false`: determines if a [`lazysort!`](/api#MolecularEvolution.lazysort!-Tuple{FelNode}) will be performed, which can reduce the amount of temporary messages that has to be initialized.
  
- `traversal=Iterators.reverse`: a function that determines the traversal, permutes an iterable.
  
- `shuffle=false`: do a randomly shuffled traversal, overrides `traversal`.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/algorithms/branchlength_optim.jl#L191-L206" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.branchlength_update!-Tuple{MolecularEvolution.UnivariateModifier, FelNode, Any}' href='#MolecularEvolution.branchlength_update!-Tuple{MolecularEvolution.UnivariateModifier, FelNode, Any}'><span class="jlbinding">MolecularEvolution.branchlength_update!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
branchlength_update!(bl_modifier::UnivariateModifier, tree::FelNode, models; <keyword arguments>)
```


A more general version of [`branchlength_optim!`](/optimization#MolecularEvolution.branchlength_optim!). Here `bl_modifier` can be either an optimizer or a sampler (or more generally, a UnivariateModifier).

**Keyword Arguments**

See [`branchlength_optim!`](/optimization#MolecularEvolution.branchlength_optim!).

::: tip Note

`bl_modifier` is a positional argument here, and not a keyword argument.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/algorithms/branchlength_optim.jl#L156-L165" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.brents_method_minimize-Tuple{Any, Real, Real, Any, Real}' href='#MolecularEvolution.brents_method_minimize-Tuple{Any, Real, Real, Any, Real}'><span class="jlbinding">MolecularEvolution.brents_method_minimize</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
brents_method_minimize(f, a::Real, b::Real, transform, t::Real; ε::Real=sqrt(eps()))
```


Brent&#39;s method for minimization.

Given a function f with a single local minimum in the interval (a,b), Brent&#39;s method returns an approximation of the x-value that minimizes f to an accuaracy between 2tol and 3tol, where tol is a combination of a relative and an absolute tolerance, tol := ε|x| + t. ε should be no smaller `2*eps`, and preferably not much less than `sqrt(eps)`, which is also the default value. eps is defined here as the machine epsilon in double precision. t should be positive.

The method combines the stability of a Golden Section Search and the superlinear convergence Successive Parabolic Interpolation has under certain conditions. The method never converges much slower than a Fibonacci search and for a sufficiently well-behaved f, convergence can be exptected to be superlinear, with an order that&#39;s usually atleast 1.3247...

**Examples**

```julia
julia> f(x) = exp(-x) - cos(x)
f (generic function with 1 method)

julia> m = brents_method_minimize(f, -1, 2, identity, 1e-7)
0.5885327257940255
```


From: Richard P. Brent, &quot;Algorithms for Minimization without Derivatives&quot; (1973). Chapter 5.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/utils/simple_optim.jl#L110-L139" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.cascading_max_state_dict-Tuple{FelNode, Any}' href='#MolecularEvolution.cascading_max_state_dict-Tuple{FelNode, Any}'><span class="jlbinding">MolecularEvolution.cascading_max_state_dict</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
cascading_max_state_dict(tree::FelNode, model; partition_list = 1:length(tree.message), node_message_dict = Dict{FelNode,Vector{<:Partition}}())
```


Takes in a tree and a model (which can be a single model, an array of models, or a function that maps FelNode-&gt;Array{&lt;:BranchModel}), and returns a dictionary mapping nodes to their inferred ancestors under the following scheme: the state that maximizes the marginal likelihood is selected at the root, and then, for each node, the maximum likelihood state is selected conditioned on the maximized state of the parent node and the observations of all descendents. A subset of partitions can be specified by partition_list, and a dictionary can be passed in to avoid re-allocating memory, in case you&#39;re running this over and over.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/algorithms/ancestors.jl#L175-L182" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.char_proportions-Tuple{Any, String}' href='#MolecularEvolution.char_proportions-Tuple{Any, String}'><span class="jlbinding">MolecularEvolution.char_proportions</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
char_proportions(seqs, alphabet::String)
```


Takes a vector of sequences and returns a vector of the proportion of each character across all sequences. An example `alphabet` argument is `MolecularEvolution.AAstring`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/discrete_models/utils/matrix_helpers.jl#L87-L92" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.collect_clades!-Tuple{FelNode, Dict{Tuple{UInt64, UInt64}, MolecularEvolution.CladeStats}, Dict{String, Int64}}' href='#MolecularEvolution.collect_clades!-Tuple{FelNode, Dict{Tuple{UInt64, UInt64}, MolecularEvolution.CladeStats}, Dict{String, Int64}}'><span class="jlbinding">MolecularEvolution.collect_clades!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Recursively collect clades from a tree.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/utils/HIPSTR.jl#L183-L185" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.collect_leaf_dists-Tuple{Vector{<:AbstractTreeNode}}' href='#MolecularEvolution.collect_leaf_dists-Tuple{Vector{<:AbstractTreeNode}}'><span class="jlbinding">MolecularEvolution.collect_leaf_dists</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
collect_leaf_dists(trees::Vector{<:AbstractTreeNode})

Returns a list of distance matrices containing the distance between the leaf nodes, which can be used to assess mixing.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/bayes/sampling.jl#L137-L141" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.colored_seq_draw-Tuple{Any, Any, AbstractString}' href='#MolecularEvolution.colored_seq_draw-Tuple{Any, Any, AbstractString}'><span class="jlbinding">MolecularEvolution.colored_seq_draw</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
colored_seq_draw(x, y, str::AbstractString; color_dict=Dict(), font_size=8pt, posx=hcenter, posy=vcenter)
```


Draw an arbitrary sequence. `color_dict` gives a mapping from characters to colors (default black). Default options for nucleotide colorings and amino acid colorings are given in the constants `NUC_COLORS` and `AA_COLORS`. This can be used along with `compose_dict` for drawing sequences at nodes in a tree (see `tree_draw`). Returns a Compose container.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/viz/tree_compose.jl#L696-L704" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.combine!-Tuple{DiscretePartition, DiscretePartition}' href='#MolecularEvolution.combine!-Tuple{DiscretePartition, DiscretePartition}'><span class="jlbinding">MolecularEvolution.combine!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
combine!(dest::P, src::P) where P<:Partition
```


Combines evidence from two partitions of the same type, storing the result in dest. Note: You should overload this for your own Partititon types.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/discrete_models/discrete_partitions.jl#L185-L190" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.copy_tree' href='#MolecularEvolution.copy_tree'><span class="jlbinding">MolecularEvolution.copy_tree</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
function copy_tree(root::FelNode, shallow_copy=false)

Returns an untangled copy of the tree. Optionally, the flag `shallow_copy` can be used to obtain a copy of the tree with only the names and branchlengths.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/nodes/FelNode.jl#L102-L107" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.count_F3x4-Tuple{Array{String}}' href='#MolecularEvolution.count_F3x4-Tuple{Array{String}}'><span class="jlbinding">MolecularEvolution.count_F3x4</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
function count_F3x4(seqs::Array{String})
```


Compute the F3x4 matrix from a set of sequences.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/discrete_models/codon_models.jl#L176-L180" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.count_F61-Tuple{Array{String}}' href='#MolecularEvolution.count_F61-Tuple{Array{String}}'><span class="jlbinding">MolecularEvolution.count_F61</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
function count_F61(seqs::Array{String}; code = universal_code)
```


Compute the F61 codon frequency vector from a set of sequences.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/discrete_models/codon_models.jl#L225-L229" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.deepequals-Union{Tuple{T}, Tuple{T, T}} where T<:AbstractTreeNode' href='#MolecularEvolution.deepequals-Union{Tuple{T}, Tuple{T, T}} where T<:AbstractTreeNode'><span class="jlbinding">MolecularEvolution.deepequals</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
deepequals(t1, t2)
```


Checks whether two trees are equal by recursively calling this on all fields, except `:parent`, in order to prevent cycles. In order to ensure that the `:parent` field is not hiding something different on both trees, ensure that each is consistent first (see: `istreeconsistent`).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/nodes/AbstractTreeNode.jl#L88-L93" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.dfs_mapreduce-Union{Tuple{T}, Tuple{AbstractTreeNode, T, Any}} where T<:Function' href='#MolecularEvolution.dfs_mapreduce-Union{Tuple{T}, Tuple{AbstractTreeNode, T, Any}} where T<:Function'><span class="jlbinding">MolecularEvolution.dfs_mapreduce</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Performs a DFS map-reduce over the tree, starting at a given node See bfs_mapreduce for more details.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/nodes/base_tree_utils.jl#L46-L49" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.discrete_name_color_dict-Tuple{AbstractTreeNode, Any}' href='#MolecularEvolution.discrete_name_color_dict-Tuple{AbstractTreeNode, Any}'><span class="jlbinding">MolecularEvolution.discrete_name_color_dict</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
discrete_name_color_dict(newt::AbstractTreeNode,tag_func; rainbow = false, scramble = false, darken = true, col_seed = nothing)
```


Takes a tree and a tag_func, which converts the leaf label into a category (ie. there should be &lt;20 of these), and returns a color dictionary that can be used to color the leaves or bubbles.

Example tag_func:     function tag_func(nam::String)         return split(nam,&quot;_&quot;)[1]     end

For prettier colors, but less discrimination: rainbow = true To randomize the rainbow color assignment: scramble = true col_seed is currently set to white, and excluded from the list of colors, to make them more visible.

Consider making your own version of this function to customize colors as you see fit.

Example use: num_leaves = 50 Ne_func(t) = 1*(e^-t).+5.0 newt = sim_tree(num_leaves,Ne_func,1.0,nstart = rand(1:num_leaves)); newt = ladderize(newt) tag_func(nam) = mod(sum(Int.(collect(nam))),7) dic = discrete_name_color_dict(newt,tag_func,rainbow = true); tree_draw(newt,line_width = 0.5mm,label_color_dict = dic)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/viz/tree_compose.jl#L511-L536" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.draw_example_tree-Tuple{}' href='#MolecularEvolution.draw_example_tree-Tuple{}'><span class="jlbinding">MolecularEvolution.draw_example_tree</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
draw_example_tree(num_leaves = 50)
```


Draws a tree and shows the code that draws it.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/viz/tree_compose.jl#L576-L580" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.endpoint_conditioned_sample_state_dict-Tuple{FelNode, Any}' href='#MolecularEvolution.endpoint_conditioned_sample_state_dict-Tuple{FelNode, Any}'><span class="jlbinding">MolecularEvolution.endpoint_conditioned_sample_state_dict</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
endpoint_conditioned_sample_state_dict(tree::FelNode, model; partition_list = 1:length(tree.message), node_message_dict = Dict{FelNode,Vector{<:Partition}}())
```


Takes in a tree and a model (which can be a single model, an array of models, or a function that maps FelNode-&gt;Array{&lt;:BranchModel}), and draws samples under the model conditions on the leaf observations. These samples are stored in the node_message_dict, which is returned. A subset of partitions can be specified by partition_list, and a dictionary can be passed in to avoid re-allocating memory, in case you&#39;re running this over and over.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/algorithms/ancestors.jl#L208-L214" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.expected_subs_per_site-Tuple{Any, Any}' href='#MolecularEvolution.expected_subs_per_site-Tuple{Any, Any}'><span class="jlbinding">MolecularEvolution.expected_subs_per_site</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
expected_subs_per_site(Q,mu)
```


Takes a rate matrix Q and an equilibrium frequency vector, and calculates the expected number of substitutions per site.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/utils/misc.jl#L144-L148" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.felsenstein!-Tuple{FelNode, Any}' href='#MolecularEvolution.felsenstein!-Tuple{FelNode, Any}'><span class="jlbinding">MolecularEvolution.felsenstein!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
felsenstein!(node::FelNode, models; partition_list = nothing)
```


Should usually be called on the root of the tree. Propagates Felsenstein pass up from the tips to the root. models can either be a single model (if the messages on the tree contain just one Partition) or an array of models, if the messages have &gt;1 Partition, or  a function that takes a node, and returns a Vector{&lt;:BranchModel} if you need the models to vary from one branch to another. partition_list (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/algorithms/felsenstein.jl#L1-L8" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.felsenstein_down!-Tuple{FelNode, Any}' href='#MolecularEvolution.felsenstein_down!-Tuple{FelNode, Any}'><span class="jlbinding">MolecularEvolution.felsenstein_down!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
felsenstein_down!(node::FelNode, models; partition_list = 1:length(tree.message), temp_message = copy_message(tree.message))
```


Should usually be called on the root of the tree. Propagates Felsenstein pass down from the root to the tips. felsenstein!() should usually be called first. models can either be a single model (if the messages on the tree contain just one Partition) or an array of models, if the messages have &gt;1 Partition, or  a function that takes a node, and returns a Vector{&lt;:BranchModel} if you need the models to vary from one branch to another. partition_list (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/algorithms/felsenstein.jl#L80-L88" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.felsenstein_roundtrip!-Tuple{FelNode, Any}' href='#MolecularEvolution.felsenstein_roundtrip!-Tuple{FelNode, Any}'><span class="jlbinding">MolecularEvolution.felsenstein_roundtrip!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
felsenstein_roundtrip!(tree::FelNode, models; partition_list = 1:length(tree.message), temp_message = copy_message(tree.message[partition_list]))
```


Should usually be called on the root of the tree. First propagates Felsenstein pass up from the tips to the root, then propagates Felsenstein pass down from the root to the tips, with the direction of time reversed (i.e. forward! = backward!). **This is useful when searching for the optimal root** (see [`root_optim!`](/optimization#MolecularEvolution.root_optim!)). models can either be a single model (if the messages on the tree contain just one Partition) or an array of models, if the messages have &gt;1 Partition, or  a function that takes a node, and returns a Vector{&lt;:BranchModel} if you need the models to vary from one branch to another. partition_list (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/algorithms/felsenstein.jl#L170-L179" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.forward!-Tuple{DiscretePartition, DiscretePartition, MolecularEvolution.PMatrixModel, FelNode}' href='#MolecularEvolution.forward!-Tuple{DiscretePartition, DiscretePartition, MolecularEvolution.PMatrixModel, FelNode}'><span class="jlbinding">MolecularEvolution.forward!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
forward!(dest::Partition, source::Partition, model::BranchModel, node::FelNode)
```


Propagate the source partition forwards along the branch to the destination partition, under the model. Note: You should overload this for your own BranchModel types.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/discrete_models/pmatrix_models/pmatrix_models.jl#L22-L27" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.gappy_Q_from_symmetric_rate_matrix-Tuple{Any, Any, Any}' href='#MolecularEvolution.gappy_Q_from_symmetric_rate_matrix-Tuple{Any, Any, Any}'><span class="jlbinding">MolecularEvolution.gappy_Q_from_symmetric_rate_matrix</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
gappy_Q_from_symmetric_rate_matrix(sym_mat, gap_rate, eq_freqs)
```


Takes a symmetric rate matrix and gap rate (governing mutations to and from gaps) and returns a gappy rate matrix. The equilibrium frequencies are multiplied on column-wise.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/discrete_models/utils/matrix_helpers.jl#L101-L106" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.get_highlighter_legend-Tuple{Any}' href='#MolecularEvolution.get_highlighter_legend-Tuple{Any}'><span class="jlbinding">MolecularEvolution.get_highlighter_legend</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
get_highlighter_legend(legend_colors)
```


Returns a Compose object given an input dictionary or pairs mapping characters to colors.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/viz/tree_compose.jl#L807-L811" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.get_max_depth-Tuple{Any, Real}' href='#MolecularEvolution.get_max_depth-Tuple{Any, Real}'><span class="jlbinding">MolecularEvolution.get_max_depth</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
get_max_depth(node,depth::Real)
```


Return the maximum depth of all children starting from the indicated node.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/viz/tree_compose.jl#L791-L795" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.get_phylo_tree-Tuple{FelNode}' href='#MolecularEvolution.get_phylo_tree-Tuple{FelNode}'><span class="jlbinding">MolecularEvolution.get_phylo_tree</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
get_phylo_tree(molev_root::FelNode; data_function = (x -> Tuple{String,Float64}[]))
```


Converts a FelNode tree to a Phylo tree. The `data_function` should return a list of tuples of the form (key, value) to be added to the Phylo tree `data` Dictionary. Any key/value pairs on the FelNode `node_data` Dict will also be added to the Phylo tree.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/viz/phylo_glue.jl#L45-L50" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.golden_section_maximize-Tuple{Any, Real, Real, Any, Real}' href='#MolecularEvolution.golden_section_maximize-Tuple{Any, Real, Real, Any, Real}'><span class="jlbinding">MolecularEvolution.golden_section_maximize</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Golden section search.

Given a function f with a single local minimum in the interval [a,b], gss returns a subset interval [c,d] that contains the minimum with d-c &lt;= tol.

**Examples**

```julia
julia> f(x) = -(x-2)^2
f (generic function with 1 method)

julia> m = golden_section_maximize(f, 1, 5, identity, 1e-10)
2.0000000000051843
```


From: https://en.wikipedia.org/wiki/Golden-section_search


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/utils/simple_optim.jl#L22-L40" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.hash_clade-Tuple{BitSet}' href='#MolecularEvolution.hash_clade-Tuple{BitSet}'><span class="jlbinding">MolecularEvolution.hash_clade</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Compute a hash for a clade based on its tips.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/utils/HIPSTR.jl#L174-L176" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.highlight_seq_draw-Tuple{Any, Any, AbstractString, Any, Any, Any}' href='#MolecularEvolution.highlight_seq_draw-Tuple{Any, Any, AbstractString, Any, Any, Any}'><span class="jlbinding">MolecularEvolution.highlight_seq_draw</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
highlight_seq_draw(x, y, str::AbstractString, region, basecolor, hicolor; fontsize=8pt, posx=hcenter, posy=vcenter)
```


Draw a sequence, highlighting the sites given in `region`. This can be used along with `compose_dict` for drawing sequences at nodes in a tree (see `tree_draw`). Returns a Compose container.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/viz/tree_compose.jl#L641-L647" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.highlighter_tree_draw-NTuple{4, Any}' href='#MolecularEvolution.highlighter_tree_draw-NTuple{4, Any}'><span class="jlbinding">MolecularEvolution.highlighter_tree_draw</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
highlighter_tree_draw(tree, ali_seqs, seqnames, master;
    highlighter_start = 1.1, highlighter_width = 1,
    coord_width = highlighter_start + highlighter_width + 0.1,
    scale_length = nothing, major_breaks = 1000, minor_breaks = 500,
    tree_args = NamedTuple[], legend_padding = 0.5cm, legend_colors = NUC_colors)
```


Draws a combined tree and highlighter plot. The vector of seqnames must match the node names in `tree`.

kwargs:
- tree_args: kwargs to pass to `tree_draw()`
  
- legend_colors: Mapping of characters to highlighter colors (default NT_colors)
  
- scale_length: Length of the scale bar
  
- highlighter_start: Canvas start for the highlighter panel
  
- highlighter_width: Canvas width for the highlighter panel
  
- coord_width: Total width of the canvas
  
- major_breaks: Numbered breaks for sequence axis
  
- minor_breaks: Ticks for sequence axis
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/viz/tree_compose.jl#L888-L907" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.internal_message_init!-Tuple{FelNode, Partition}' href='#MolecularEvolution.internal_message_init!-Tuple{FelNode, Partition}'><span class="jlbinding">MolecularEvolution.internal_message_init!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
internal_message_init!(tree::FelNode, partition::Partition)

Initializes the message template for each node in the tree, as an array of the partition.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/nodes/FelNode.jl#L70-L74" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.internal_message_init!-Tuple{FelNode, Vector{<:Partition}}' href='#MolecularEvolution.internal_message_init!-Tuple{FelNode, Vector{<:Partition}}'><span class="jlbinding">MolecularEvolution.internal_message_init!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
internal_message_init!(tree::FelNode, empty_message::Vector{<:Partition})

Initializes the message template for each node in the tree, allocating space for each partition.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/nodes/FelNode.jl#L55-L59" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.internal_nodes-Tuple{Any}' href='#MolecularEvolution.internal_nodes-Tuple{Any}'><span class="jlbinding">MolecularEvolution.internal_nodes</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
internal_nodes(tree)
```


Returns the internal nodes of the tree (including the root), as a vector.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/simple_interface.jl#L60-L64" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.istreeconsistent-Tuple{T} where T<:AbstractTreeNode' href='#MolecularEvolution.istreeconsistent-Tuple{T} where T<:AbstractTreeNode'><span class="jlbinding">MolecularEvolution.istreeconsistent</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
istreeconsistent(root)
```


Checks whether the `:parent` field is set to be consistent with the `:child` field for all nodes in the subtree. 


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/nodes/AbstractTreeNode.jl#L61-L65" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.ladder_tree_sim-Tuple{Any}' href='#MolecularEvolution.ladder_tree_sim-Tuple{Any}'><span class="jlbinding">MolecularEvolution.ladder_tree_sim</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
ladder_tree_sim(ntaxa)
```


Simulates a ladder-like tree, using constant population size but heterochronous sampling, under a coalescent model.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/sim_tree.jl#L154-L158" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.lazyprep!-Tuple{FelNode, Vector{<:Partition}}' href='#MolecularEvolution.lazyprep!-Tuple{FelNode, Vector{<:Partition}}'><span class="jlbinding">MolecularEvolution.lazyprep!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
lazyprep!(tree::FelNode, initial_message::Vector{<:Partition}; partition_list = 1:length(tree.message), direction::LazyDirection = LazyUp())
```


Extra, intermediate step of tree preparations between initializing messages across the tree and calling message passing algorithms with `LazyPartition`.
1. Perform a `lazysort!` on tree to obtain the optimal tree for a lazy `felsenstein!` prop, or a `sample_down!`.
  
2. Fix `tree.parent_message` to an initial message.
  
3. Preallocate sufficiently many inner partitions needed for a `felsenstein!` prop, or a `sample_down!`.
  
4. Specialized preparations based on the direction of the operations (`forward!`, `backward!`). `LazyDown` or `LazyUp`.
  

See also `LazyDown`, `LazyUp`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/lazy_models/lazy_partition.jl#L200-L210" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.lazysort!-Tuple{FelNode}' href='#MolecularEvolution.lazysort!-Tuple{FelNode}'><span class="jlbinding">MolecularEvolution.lazysort!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>


- Should be run on a tree containing LazyPartitions before running `felsenstein!`. Sorts for a minimal count of active partitions during a felsenstein!
  
- Returns the minimum length of memoryblocks (-1) required for a `felsenstein!` prop. We need a temporary memoryblock during `backward!`, hence the &#39;-1&#39;.
  

::: tip Note

Since felsenstein! uses a stack, we want to avoid having long node.children[1].children[1]... chains

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/lazy_models/lazy_partition.jl#L152-L157" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.leaf_distmat-Tuple{Any}' href='#MolecularEvolution.leaf_distmat-Tuple{Any}'><span class="jlbinding">MolecularEvolution.leaf_distmat</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
leaf_distmat(tree)
```


Returns a matrix of the distances between the leaf nodes where the index on the columns and rows are sorted by the leaf names.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/bayes/sampling.jl#L150-L154" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.leaf_names-Tuple{FelNode}' href='#MolecularEvolution.leaf_names-Tuple{FelNode}'><span class="jlbinding">MolecularEvolution.leaf_names</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
leaf_names(tree)
```


Returns the names of the leaves of the tree.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/simple_interface.jl#L24-L28" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.leaf_samples-Tuple{FelNode}' href='#MolecularEvolution.leaf_samples-Tuple{FelNode}'><span class="jlbinding">MolecularEvolution.leaf_samples</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
leaf_samples(tree; partition_inds = 1)
```


Returns the result of `partition2obs` for each leaf of the tree. Can be used eg. after `sample_down!` is called. If using a eg. codon model, this will extract a string from the CodonPartition on each leaf. Acts upon the first partition by default, but this can be changed by setting `partition_inds`, which can also be a vector of indices,  in which case the result will be a vector for each leaf.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/simple_interface.jl#L4-L11" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.leaves-Tuple{Any}' href='#MolecularEvolution.leaves-Tuple{Any}'><span class="jlbinding">MolecularEvolution.leaves</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
leaves(tree)
```


Returns the leaves of the tree, as a vector.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/simple_interface.jl#L46-L50" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.linear_scale-NTuple{5, Any}' href='#MolecularEvolution.linear_scale-NTuple{5, Any}'><span class="jlbinding">MolecularEvolution.linear_scale</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
linear_scale(val,in_min,in_max,out_min,out_max)
```


Linearly maps val which lives in [in_min,in_max] to a value in [out_min,out_max]


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/utils/misc.jl#L199-L203" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.log_likelihood!-Tuple{FelNode, Any}' href='#MolecularEvolution.log_likelihood!-Tuple{FelNode, Any}'><span class="jlbinding">MolecularEvolution.log_likelihood!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
log_likelihood!(tree::FelNode, models; partition_list = nothing)
```


First re-computes the upward felsenstein pass, and then computes the log likelihood of this tree. models can either be a single model (if the messages on the tree contain just one Partition) or an array of models, if the messages have &gt;1 Partition, or  a function that takes a node, and returns a Vector{&lt;:BranchModel} if you need the models to vary from one branch to another. partition_list (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/algorithms/lls.jl#L45-L52" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.log_likelihood-Tuple{FelNode, BranchModel}' href='#MolecularEvolution.log_likelihood-Tuple{FelNode, BranchModel}'><span class="jlbinding">MolecularEvolution.log_likelihood</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
log_likelihood(tree::FelNode, models; partition_list = nothing)
```


Computed the log likelihood of this tree. Requires felsenstein!() to have been run. models can either be a single model (if the messages on the tree contain just one Partition) or an array of models, if the messages have &gt;1 Partition, or  a function that takes a node, and returns a Vector{&lt;:BranchModel} if you need the models to vary from one branch to another. partition_list (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/algorithms/lls.jl#L27-L34" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.longest_path-Tuple{FelNode}' href='#MolecularEvolution.longest_path-Tuple{FelNode}'><span class="jlbinding">MolecularEvolution.longest_path</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Returns the longest path in a tree For convenience, this is returned as two lists of form:     [leaf_node, parent_node, .... root] Where the leaf_node nodes are selected to be the furthest away


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/nodes/base_tree_utils.jl#L106-L111" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.marginal_state_dict-Tuple{FelNode, Any}' href='#MolecularEvolution.marginal_state_dict-Tuple{FelNode, Any}'><span class="jlbinding">MolecularEvolution.marginal_state_dict</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
marginal_state_dict(tree::FelNode, model; partition_list = 1:length(tree.message), node_message_dict = Dict{FelNode,Vector{<:Partition}}())
```


Takes in a tree and a model (which can be a single model, an array of models, or a function that maps FelNode-&gt;Array{&lt;:BranchModel}), and returns a dictionary mapping nodes to their marginal reconstructions (ie. P(state|all observations,model)). A subset of partitions can be specified by partition_list, and a dictionary can be passed in to avoid re-allocating memory, in case you&#39;re running this over and over.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/algorithms/ancestors.jl#L111-L117" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.matrix_for_display-Tuple{Any, Any}' href='#MolecularEvolution.matrix_for_display-Tuple{Any, Any}'><span class="jlbinding">MolecularEvolution.matrix_for_display</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
matrix_for_display(Q,labels)
```


Takes a numerical matrix and a vector of labels, and returns a typically mixed type matrix with the numerical values and the labels. This is to easily visualize rate matrices in eg. the REPL.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/discrete_models/utils/matrix_helpers.jl#L196-L201" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.metropolis_sample-Tuple{AbstractUpdate, FelNode, Any, Any}' href='#MolecularEvolution.metropolis_sample-Tuple{AbstractUpdate, FelNode, Any, Any}'><span class="jlbinding">MolecularEvolution.metropolis_sample</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
metropolis_sample(
    update!::AbstractUpdate,
    initial_tree::FelNode,
    models,
    num_of_samples;
    partition_list = 1:length(initial_tree.message),
    burn_in = 1000,
    sample_interval = 10,
    collect_LLs = false,
    collect_models = false,
    midpoint_rooting = false,
    ladderize = false,
)
```


Samples tree topologies from a posterior distribution using a custom `update!` function.

**Arguments**
- `update!`: A callable that takes (tree::FelNode, models) and updates `tree` and `models`. One call to `update!` corresponds to one iteration of the Metropolis algorithm.
  
- `initial_tree`: An initial tree topology with the leaves populated with data, for the likelihood calculation.
  
- `models`: A list of branch models.
  
- `num_of_samples`: The number of tree samples drawn from the posterior.
  
- `partition_list`: (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over (but you probably want to sample with all partitions, the default option).
  
- `burn_in`: The number of samples discarded at the start of the Markov Chain.
  
- `sample_interval`: The distance between samples in the underlying Markov Chain (to reduce sample correlation).
  
- `collect_LLs`: Specifies if the function should return the log-likelihoods of the trees.
  
- `collect_models`: Specifies if the function should return the models.
  
- `midpoint_rooting`: Specifies whether the drawn samples should be midpoint rerooted (Important! Should only be used for time-reversible branch models starting in equilibrium).
  

::: tip Note

The leaves of the initial tree should be populated with data and felsenstein! should be called on the initial tree before calling this function.

:::

**Returns**
- `samples`: The trees drawn from the posterior. Returns shallow tree copies, which needs to be repopulated before running felsenstein! etc. 
  
- `sample_LLs`: The associated log-likelihoods of the tree (optional).
  
- `sample_models`: The models drawn from the posterior (optional). The models can be collapsed into it&#39;s parameters with `collapse_models`.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/bayes/sampling.jl#L1-L37" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.metropolis_sample-Tuple{FelNode, Vector{<:BranchModel}, Any}' href='#MolecularEvolution.metropolis_sample-Tuple{FelNode, Vector{<:BranchModel}, Any}'><span class="jlbinding">MolecularEvolution.metropolis_sample</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
metropolis_sample(
    initial_tree::FelNode,
    models::Vector{<:BranchModel},
    num_of_samples;
    bl_sampler::UnivariateSampler = BranchlengthSampler(Normal(0,2), Normal(-1,1))
    burn_in=1000, 
    sample_interval=10,
    collect_LLs = false,
    midpoint_rooting=false,
)
```


A convenience method. One step of the Metropolis algorithm is performed by calling [`nni_update!`](/api#MolecularEvolution.nni_update!-Tuple{Function,%20FelNode,%20Any}) with `softmax_sampler` and [`branchlength_update!`](/api#MolecularEvolution.branchlength_update!-Tuple{MolecularEvolution.UnivariateModifier,%20FelNode,%20Any}) with `bl_sampler`.

**Additional Arguments**
- `bl_sampler`: Sampler used to drawn branchlengths from the posterior.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/bayes/sampling.jl#L108-L124" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.metropolis_step-Tuple{Function, Any, Any}' href='#MolecularEvolution.metropolis_step-Tuple{Function, Any, Any}'><span class="jlbinding">MolecularEvolution.metropolis_step</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
metropolis_step(LL::Function, modifier, curr_value)
```


Does a standard metropolis step in an MCMC, i.e. proposes a candidate symmetrically and returns the next state in the chain, decided by the candidate being rejected or not.

**Interface**

You need a `MySampler <: Any` to implement
- `proposal(modifier::MySampler, curr_value)`
  
- `log_prior(modifier::MySampler, x)`
  
- `apply_decision(modifier::MySampler, accept::Bool)`
  

`LL` is by default called on `curr_value` and the returned value of `proposal`. Although, it is possible to transform the current value before proposing a new value, and then take the inverse transform to match the argument `LL` expects.

**Extended interface**

**Hastings**

To allow for asymmetric proposals, you must overload
- `log_proposal(modifier::MySampler, x, conditioned_on)`
  

which returns a constant (`0.0` in particular) by default.

**Transformations**

To make proposals in a transformed space, you overload
- `tr(modifier::MySampler, x)`
  
- `invtr(modifier::MySampler, x)`
  

which are identity transformations by default.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/utils/simple_sample.jl#L27-L51" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.midpoint-Tuple{FelNode}' href='#MolecularEvolution.midpoint-Tuple{FelNode}'><span class="jlbinding">MolecularEvolution.midpoint</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Returns a midpoint as a node and a distance above it where the midpoint is


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/nodes/base_tree_utils.jl#L126" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.mix-Union{Tuple{SWMPartition{PType}}, Tuple{PType}} where PType<:DiscretePartition' href='#MolecularEvolution.mix-Union{Tuple{SWMPartition{PType}}, Tuple{PType}} where PType<:DiscretePartition'><span class="jlbinding">MolecularEvolution.mix</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
mix(swm_part::SWMPartition{PType} ) where {PType <: MultiSitePartition}
```


`mix` collapses a Site-Wise Mixture partition to a single component partition, weighted by the site-wise likelihoods for each component, and the init weights. Specifically, it takes a `SWMPartition{Ptype}` and returns a `PType`. You&#39;ll need to have this implemented for certain helper functionality if you&#39;re playing with new kinds of SWMPartitions that aren&#39;t mixtures of `DiscretePartitions`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/compound_models/swm.jl#L188-L194" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.name2node_dict-Tuple{AbstractTreeNode}' href='#MolecularEvolution.name2node_dict-Tuple{AbstractTreeNode}'><span class="jlbinding">MolecularEvolution.name2node_dict</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



name2node_dict(root)

Returns a dictionary of leaf nodes, indexed by node.name. Can be used to associate sequences with leaf nodes.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/nodes/AbstractTreeNode.jl#L540-L544" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.newick-Tuple{AbstractTreeNode}' href='#MolecularEvolution.newick-Tuple{AbstractTreeNode}'><span class="jlbinding">MolecularEvolution.newick</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
newick(root)
```


Returns a newick string representation of the tree.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/nodes/AbstractTreeNode.jl#L626-L630" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.nni_optim!-Tuple' href='#MolecularEvolution.nni_optim!-Tuple'><span class="jlbinding">MolecularEvolution.nni_optim!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
nni_optim!(tree::FelNode, models; <keyword arguments>)
```


Considers local branch swaps for all branches recursively, maintaining the integrity of the messages. Requires felsenstein!() to have been run first. models can either be a single model (if the messages on the tree contain just one Partition) or an array of models, if the messages have &gt;1 Partition, or  a function that takes a node, and returns a Vector{&lt;:BranchModel} if you need the models to vary from one branch to another.

**Keyword Arguments**
- `partition_list=nothing`: (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over (but you probably want to optimize tree topology with all models, the default option).
  
- `selection_rule = x -> argmax(x)`: a function that takes the current and proposed log likelihoods and selects a nni configuration. Note that the current log likelihood is stored at x[1].
  
- `sort_tree=false`: determines if a [`lazysort!`](/api#MolecularEvolution.lazysort!-Tuple{FelNode}) will be performed, which can reduce the amount of temporary messages that has to be initialized.
  
- `traversal=Iterators.reverse`: a function that determines the traversal, permutes an iterable.
  
- `shuffle=false`: do a randomly shuffled traversal, overrides `traversal`.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/algorithms/nni_optim.jl#L323-L337" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.nni_update!-Tuple{Function, FelNode, Any}' href='#MolecularEvolution.nni_update!-Tuple{Function, FelNode, Any}'><span class="jlbinding">MolecularEvolution.nni_update!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
nni_update!(selection_rule::Function, tree::FelNode, models; <keyword arguments>)
```


A more verbose version of [`nni_optim!`](/optimization#MolecularEvolution.nni_optim!).

**Keyword Arguments**

See [`nni_optim!`](/optimization#MolecularEvolution.nni_optim!).

::: tip Note

`selection_rule` is a positional argument here, and not a keyword argument.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/algorithms/nni_optim.jl#L287-L296" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.node_distances-Tuple{AbstractTreeNode}' href='#MolecularEvolution.node_distances-Tuple{AbstractTreeNode}'><span class="jlbinding">MolecularEvolution.node_distances</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Compute the distance to all other nodes from a given node


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/nodes/base_tree_utils.jl#L58" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.node_names-Tuple{FelNode}' href='#MolecularEvolution.node_names-Tuple{FelNode}'><span class="jlbinding">MolecularEvolution.node_names</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
node_names(tree)
```


Returns the names of the nodes of the tree.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/simple_interface.jl#L31-L35" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.node_samples-Tuple{FelNode}' href='#MolecularEvolution.node_samples-Tuple{FelNode}'><span class="jlbinding">MolecularEvolution.node_samples</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
node_samples(tree; partition_inds = 1)
```


Returns the result of `partition2obs` for each node of the tree (including internal nodes, and the root). Can be used eg. after `sample_down!` is called. If using a eg. codon model, this will extract a string from the CodonPartition on each node. Acts upon the first partition by default, but this can be changed by setting `partition_inds`, which can also be a vector of indices,  in which case the result will be a vector for each node.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/simple_interface.jl#L14-L21" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.nodes-Tuple{Any}' href='#MolecularEvolution.nodes-Tuple{Any}'><span class="jlbinding">MolecularEvolution.nodes</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
nodes(tree)
```


Returns the nodes of the tree (including internal nodes and the root), as a vector.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/simple_interface.jl#L53-L57" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.nonreversibleQ-Tuple{Any}' href='#MolecularEvolution.nonreversibleQ-Tuple{Any}'><span class="jlbinding">MolecularEvolution.nonreversibleQ</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
nonreversibleQ(param_vec)
```


Takes a vector of parameters and returns a nonreversible rate matrix.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/discrete_models/utils/matrix_helpers.jl#L151-L155" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.parent_list-Tuple{FelNode}' href='#MolecularEvolution.parent_list-Tuple{FelNode}'><span class="jlbinding">MolecularEvolution.parent_list</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Provides a list of parent nodes nodes from this node up to the root node


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/nodes/base_tree_utils.jl#L73" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.partition2obs-Tuple{DiscretePartition, String}' href='#MolecularEvolution.partition2obs-Tuple{DiscretePartition, String}'><span class="jlbinding">MolecularEvolution.partition2obs</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
partition2obs(part::Partition)
```


Extracts the most likely state from a Partition, transforming it into a convenient type. For example, a NucleotidePartition will be transformed into a nucleotide sequence of type String. Note: You should overload this for your own Partititon types.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/discrete_models/utils/seq_to_vec.jl#L70-L76" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.plot_multiple_trees-Tuple{Any, Any}' href='#MolecularEvolution.plot_multiple_trees-Tuple{Any, Any}'><span class="jlbinding">MolecularEvolution.plot_multiple_trees</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plot_multiple_trees(trees, inf_tree; <keyword arguments>)
```


Plots multiple phylogenetic trees against a reference tree, `inf_tree`. For each **tree** in `trees`, a linear Weighted Least Squares (WLS) problem (parameterized by the `weight_fn` keyword) is solved for the x-positions of the matching nodes between `inf_tree` and **tree**.

**Keyword Arguments**
- `node_size=4`: the size of the nodes in the plot.
  
- `line_width=0.5`: the width of the branches from `trees`.
  
- `font_size=10`: the font size for the leaf labels.
  
- `margin=1.5`: the margin between a leaf node and its label.
  
- `line_alpha=0.05`: the transparency level of the branches from `trees`.
  
- `y_jitter=0.0`: the standard deviation of the noise in the y-coordinate.
  
- `weight_fn=n::FelNode -> ifelse(isroot(n), 1.0, 0.0))`: a function that assigns a weight to a node for the WLS problem.
  
- `opt_scale=true`: whether to include a scaling parameter for the WLS problem.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/viz/multiple_trees.jl#L74-L89" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.populate_tree!-Tuple{FelNode, Partition, Any, Any}' href='#MolecularEvolution.populate_tree!-Tuple{FelNode, Partition, Any, Any}'><span class="jlbinding">MolecularEvolution.populate_tree!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
populate_tree!(tree::FelNode, starting_message, names, data; init_all_messages = true, tolerate_missing = 1, leaf_name_transform = x -> x)
```


Takes a tree, and a `starting_message` (which will serve as the memory template for populating messages all over the tree). `starting_message` can be a message (ie. a vector of Partitions), but will also work with a single Partition (although the tree) will still be populated with a length-1 vector of Partitions. Further, as long as `obs2partition` is implemented for your Partition type, the leaf nodes will be populated with the data from `data`, matching the names on each leaf. When a leaf on the tree has a name that doesn&#39;t match anything in `names`, then if
- `tolerate_missing = 0`, an error will be thrown
  
- `tolerate_missing = 1`, a warning will be thrown, and the message will be set to the uninformative message (requires identity!(::Partition) to be defined)
  
- `tolerate_missing = 2`, the message will be set to the uninformative message, without warnings (requires identity!(::Partition) to be defined)
  

A renaming function that can eg. strip tags from the tree when matching leaf names with `names` can be passed to `leaf_name_transform`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/utils/misc.jl#L110-L122" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.promote_internal-Tuple{FelNode}' href='#MolecularEvolution.promote_internal-Tuple{FelNode}'><span class="jlbinding">MolecularEvolution.promote_internal</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
promote_internal(tree::FelNode)
```


Creates a new tree similar to the given tree, but with &#39;dummy&#39; leaf nodes (w/ zero branchlength) representing each internal node (for drawing / evenly spacing labels internal nodes).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/viz/tree_compose.jl#L615-L620" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.quadratic_CI-Tuple{Function, Vector, Int64}' href='#MolecularEvolution.quadratic_CI-Tuple{Function, Vector, Int64}'><span class="jlbinding">MolecularEvolution.quadratic_CI</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
quadratic_CI(f::Function,opt_params::Vector, param_ind::Int; rate_conf_level = 0.99, nudge_amount = 0.01)
```


Takes a NEGATIVE log likelihood function (compatible with Optim.jl), a vector of maximizing parameters, an a parameter index. Returns the quadratic confidence interval.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/utils/simple_optim.jl#L357-L362" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.quadratic_CI-Tuple{Vector, Vector}' href='#MolecularEvolution.quadratic_CI-Tuple{Vector, Vector}'><span class="jlbinding">MolecularEvolution.quadratic_CI</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
quadratic_CI(xvec,yvec; rate_conf_level = 0.99)
```


Takes xvec, a vector of parameter values, and yvec, a vector of log likelihood evaluations (note: NOT the negative LLs you) might use with Optim.jl. Returns the confidence intervals computed by a quadratic approximation to the LL.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/utils/simple_optim.jl#L333-L339" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.read_fasta-Tuple{String}' href='#MolecularEvolution.read_fasta-Tuple{String}'><span class="jlbinding">MolecularEvolution.read_fasta</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
read_fasta(filepath::String)
```


Reads in a fasta file and returns a tuple of (seqnames, seqs).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/utils/fasta_io.jl#L4-L8" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.read_newick_tree-Tuple{String}' href='#MolecularEvolution.read_newick_tree-Tuple{String}'><span class="jlbinding">MolecularEvolution.read_newick_tree</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



read_newick_tree(treefile)

Reads in a tree from a file, of type FelNode


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/utils/misc.jl#L256-L260" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.refresh!-Tuple{FelNode, Any}' href='#MolecularEvolution.refresh!-Tuple{FelNode, Any}'><span class="jlbinding">MolecularEvolution.refresh!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
refresh!(tree::FelNode, models; partition_list = 1:length(tree.message))
```


Run `felsenstein!` and `felsenstein_down!` on `tree` with `models` to refresh messages.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/algorithms/AbstractUpdate.jl#L144-L148" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.reversibleQ-Tuple{Any, Any}' href='#MolecularEvolution.reversibleQ-Tuple{Any, Any}'><span class="jlbinding">MolecularEvolution.reversibleQ</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
reversibleQ(param_vec,eq_freqs)
```


Takes a vector of parameters and equilibrium frequencies and returns a reversible rate matrix. The parameters are the upper triangle of the rate matrix, with the diagonal elements omitted, and the equilibrium frequencies are multiplied column-wise.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/discrete_models/utils/matrix_helpers.jl#L122-L128" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.root2tip_distances-Tuple{AbstractTreeNode}' href='#MolecularEvolution.root2tip_distances-Tuple{AbstractTreeNode}'><span class="jlbinding">MolecularEvolution.root2tip_distances</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
root2tips(root::AbstractTreeNode)
```


Returns a vector of root-to-tip distances, and a node-to-index dictionary. Be aware that this dictionary will break when any of the node content (ie. anything on the tree) changes.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/nodes/AbstractTreeNode.jl#L773-L778" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.root_optim!-Tuple{FelNode, Any}' href='#MolecularEvolution.root_optim!-Tuple{FelNode, Any}'><span class="jlbinding">MolecularEvolution.root_optim!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
root_optim!(tree::FelNode, models; <keyword arguments>)
```


Optimizes the root position and root state of a tree. Returns the new, optimal root node. models can either be a single model (if the messages on the tree contain just one Partition) or an array of models, if the messages have &gt;1 Partition, or  a function that takes a node, and returns a Vector{&lt;:BranchModel} if you need the models to vary from one branch to another.

**Keyword Arguments**
- `partition_list=1:length(tree.message)`: (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over (but you probably want to optimize root position and root state with all models, the default option).
  
- `root_LL!=default_root_LL_wrapper(tree.parent_message[partition_list])`: a function that takes a message and returns a (optimal) parent message and LL (log likelihood). The default option uses the constant `tree.parent_message[partition_list]` as parent message for all root-candidates.
  
- `K=10`: the number of equidistant root-candidate points along a branch. (only to be used in the frequentist framework!?)
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/algorithms/root_optim.jl#L108-L119" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.root_update!-Tuple{RootUpdate, FelNode, Any}' href='#MolecularEvolution.root_update!-Tuple{RootUpdate, FelNode, Any}'><span class="jlbinding">MolecularEvolution.root_update!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
root_update!(root_update::RootUpdate, tree::FelNode, models; partition_list = 1:length(tree.message))
```


A more general version of [`root_optim!`](/optimization#MolecularEvolution.root_optim!). Here `root_update` can be either an optimization or a sampling (or more generally, a RootUpdate).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/algorithms/root_optim.jl#L88-L92" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.sample_down!-Tuple{FelNode, Any, Any}' href='#MolecularEvolution.sample_down!-Tuple{FelNode, Any, Any}'><span class="jlbinding">MolecularEvolution.sample_down!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



sample_down!(root::FelNode,models,partition_list)

Generates samples under the model. The root.parent_message is taken as the starting distribution, and node.message contains the sampled messages. models can either be a single model (if the messages on the tree contain just one Partition) or an array of models, if the messages have &gt;1 Partition, or  a function that takes a node, and returns a Vector{&lt;:BranchModel} if you need the models to vary from one branch to another. partition_list (eg. 1:3 or [1,3,5]) lets you choose which partitions to run over.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/algorithms/generative.jl#L14-L21" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.sample_from_message!-Tuple{Vector{<:Partition}}' href='#MolecularEvolution.sample_from_message!-Tuple{Vector{<:Partition}}'><span class="jlbinding">MolecularEvolution.sample_from_message!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
sample_from_message!(message::Vector{<:Partition})
```


#Replaces an uncertain message with a sample from the distribution represented by each partition.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/algorithms/generative.jl#L5-L9" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.savefig_tweakSVG-Tuple{Any, Context}' href='#MolecularEvolution.savefig_tweakSVG-Tuple{Any, Context}'><span class="jlbinding">MolecularEvolution.savefig_tweakSVG</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
savefig_tweakSVG(fname, plot::Context; width = 10cm, height = 10cm, linecap_round = true, white_background = true)
```


Saves a figure created using the `Compose` approach, but tweaks the SVG after export.

eg. `savefig_tweakSVG("export.svg",pl)`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/viz/tree_compose.jl#L1012-L1018" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.savefig_tweakSVG-Tuple{Any, Plots.Plot}' href='#MolecularEvolution.savefig_tweakSVG-Tuple{Any, Plots.Plot}'><span class="jlbinding">MolecularEvolution.savefig_tweakSVG</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
savefig_tweakSVG(fname, plot::Plots.Plot; hack_bounding_box = true, new_viewbox = nothing, linecap_round = true)
```


Note: Might only work if you&#39;re using the GR backend!! Saves a figure created using the `Phylo` `Plots` recipe, but tweaks the SVG after export. `new_viewbox` needs to be an array of 4 numbers, typically starting at `[0 0 plot_width*4 plot_height*4]` but this lets you add shifts, in case the plot is getting cut off.

eg. `savefig_tweakSVG("export.svg",pl, new_viewbox = [-100, -100, 3000, 4500])`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/viz/phylo_glue.jl#L84-L93" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.shortest_path_between_nodes-Tuple{FelNode, FelNode}' href='#MolecularEvolution.shortest_path_between_nodes-Tuple{FelNode, FelNode}'><span class="jlbinding">MolecularEvolution.shortest_path_between_nodes</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Shortest path between nodes, returned as two lists, each starting with one of the two nodes,  and ending with the common ancestor


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/nodes/base_tree_utils.jl#L84-L87" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.sibling_inds-Tuple{AbstractTreeNode}' href='#MolecularEvolution.sibling_inds-Tuple{AbstractTreeNode}'><span class="jlbinding">MolecularEvolution.sibling_inds</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



sibling_inds(node)

Returns logical indices of the siblings in the parent&#39;s child&#39;s vector.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/nodes/AbstractTreeNode.jl#L526-L530" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.siblings-Tuple{AbstractTreeNode}' href='#MolecularEvolution.siblings-Tuple{AbstractTreeNode}'><span class="jlbinding">MolecularEvolution.siblings</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



siblings(node)

Returns a vector of siblings of node.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/nodes/AbstractTreeNode.jl#L512-L516" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.sim_tree-Tuple{Int64, Any, Any}' href='#MolecularEvolution.sim_tree-Tuple{Int64, Any, Any}'><span class="jlbinding">MolecularEvolution.sim_tree</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
sim_tree(add_limit::Int,Ne_func,sample_rate_func; nstart = 1, time = 0.0, mutation_rate = 1.0, T = Float64)
```


Simulates a tree of type FelNode{T}. Allows an effective population size function (Ne_func), as well as a sample rate function (sample_rate_func), which can also just be constants.

Ne_func(t) = (sin(t/10)+1)*100.0 + 10.0 root = sim_tree(600,Ne_func,1.0) simple_tree_draw(ladderize(root))


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/sim_tree.jl#L94-L103" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.sim_tree-Tuple{}' href='#MolecularEvolution.sim_tree-Tuple{}'><span class="jlbinding">MolecularEvolution.sim_tree</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
sim_tree(;n = 10)
```


Simulates tree with constant population size.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/sim_tree.jl#L132-L137" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.simple_radial_tree_plot-Tuple{FelNode}' href='#MolecularEvolution.simple_radial_tree_plot-Tuple{FelNode}'><span class="jlbinding">MolecularEvolution.simple_radial_tree_plot</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
simple_radial_tree_plot(root::FelNode; canvas_width = 10cm, line_color = "black", line_width = 0.1mm)
```


Draws a radial tree. No frills. No labels. Canvas height is automatically determined to avoid distorting the tree.

newt = better_newick_import(&quot;((A:1,B:1,C:1,D:1,E:1,F:1,G:1):1,(H:1,I:1):1);&quot;, FelNode{Float64}); simple_radial_tree_plot(newt,line_width = 0.5mm,root_angle = 7/10)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/viz/tree_compose.jl#L423-L430" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.simple_tree_draw-Tuple{FelNode}' href='#MolecularEvolution.simple_tree_draw-Tuple{FelNode}'><span class="jlbinding">MolecularEvolution.simple_tree_draw</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



img = simple_tree_draw(tree::FelNode; canvas_width = 15cm, canvas_height = 15cm, line_color = &quot;black&quot;, line_width = 0.1mm)

A line drawing of a tree with very few options.

```julia
img = simple_tree_draw(tree)
img |> SVG("imgout.svg",10cm, 10cm)
OR
using Cairo
img |> PDF("imgout.pdf",10cm, 10cm)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/viz/tree_compose.jl#L92-L104" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.standard_tree_sim-Tuple{Any}' href='#MolecularEvolution.standard_tree_sim-Tuple{Any}'><span class="jlbinding">MolecularEvolution.standard_tree_sim</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
standard_tree_sim(ntaxa)
```


Simulates a tree with logistic population growth, under a coalescent model.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/sim_tree.jl#L143-L147" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.total_LL-Tuple{Partition}' href='#MolecularEvolution.total_LL-Tuple{Partition}'><span class="jlbinding">MolecularEvolution.total_LL</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



total_LL(p::Partition)

If called on the root, it returns the log likelihood associated with that partition. Can be overloaded for complex partitions without straightforward site log likelihoods.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/algorithms/lls.jl#L2-L7" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.tree2distances-Tuple{AbstractTreeNode}' href='#MolecularEvolution.tree2distances-Tuple{AbstractTreeNode}'><span class="jlbinding">MolecularEvolution.tree2distances</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
tree2distances(root::AbstractTreeNode)
```


Returns a distance matrix for all pairs of leaf nodes, and a node-to-index dictionary. Be aware that this dictionary will break when any of the node content (ie. anything on the tree) changes.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/nodes/AbstractTreeNode.jl#L736-L741" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.tree2shared_branch_lengths-Tuple{AbstractTreeNode}' href='#MolecularEvolution.tree2shared_branch_lengths-Tuple{AbstractTreeNode}'><span class="jlbinding">MolecularEvolution.tree2shared_branch_lengths</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
tree2distances(root::AbstractTreeNode)
```


Returns a distance matrix for all pairs of leaf nodes, and a node-to-index dictionary. Be aware that this dictionary will break when any of the node content (ie. anything on the tree) changes.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/nodes/AbstractTreeNode.jl#L791-L796" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.tree_draw-Tuple{FelNode}' href='#MolecularEvolution.tree_draw-Tuple{FelNode}'><span class="jlbinding">MolecularEvolution.tree_draw</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
tree_draw(tree::FelNode;
    canvas_width = 15cm, canvas_height = 15cm,
    stretch_for_labels = 2.0, draw_labels = true,
    line_width = 0.1mm, font_size = 4pt,
    min_dot_size = 0.00, max_dot_size = 0.01,
    line_opacity = 1.0,
    dot_opacity = 1.0,
    name_opacity = 1.0,
    horizontal = true,
    dot_size_dict = Dict(), dot_size_default = 0.0,
    dot_color_dict = Dict(), dot_color_default = "black",
    line_color_dict = Dict(), line_color_default = "black",
    label_color_dict = Dict(), label_color_default = "black",
    nodelabel_dict = Dict(),compose_dict = Dict()
    )
```


Draws a tree with a number of self-explanatory options. Dictionaries that map a node to a color/size are used to control per-node plotting options. `compose_dict` must be a `FelNode->function(x,y)` dictionary that returns a `compose()` struct.

Example using `compose_dict`

```julia
str_tree = "(((((tax24:0.09731668728575642,(tax22:0.08792233964843627,tax18:0.9210388482867483):0.3200367900275155):0.6948314526087965,(tax13:1.9977212308725611,(tax15:0.4290074347886068,(tax17:0.32928401808187824,(tax12:0.3860215462534818,tax16:0.2197134841232339):0.1399122681886174):0.05744611946245004):1.4686085778061146):0.20724159879522402):0.4539334554156126,tax28:0.4885576926440158):0.002162260013924424,tax26:0.9451873777301325):3.8695419798779387,((tax29:0.10062813251515536,tax27:0.27653633028085006):0.04262434258357507,(tax25:0.009345653929737636,((tax23:0.015832941547076644,(tax20:0.5550597590956172,((tax8:0.6649025646927402,tax9:0.358506423199849):0.1439516404012261,tax11:0.01995439013213013):1.155181296134081):0.17930021667907567):0.10906638146207207,((((((tax6:0.013708993438720255,tax5:0.061144001556547097):0.1395453591567641,tax3:0.4713722705245479):0.07432598428904214,tax1:0.5993347898257291):1.0588025698844894,(tax10:0.13109032492533992,(tax4:0.8517302241963356,(tax2:0.8481963081549965,tax7:0.23754095940676642):0.2394313086297733):0.43596704123297675):0.08774657269409454):0.9345533723114966,(tax14:0.7089558245245173,tax19:0.444897137240675):0.08657675809803095):0.01632062723968511,tax21:0.029535281963725537):0.49502691718938285):0.25829576024240986):0.7339777396780424):4.148878039524972):0.0"
newt = gettreefromnewick(str_tree, FelNode)
ladderize!(newt)
compose_dict = Dict()
for n in getleaflist(newt)
    #Replace the rand(4) with the frequencies you actually want.
    compose_dict[n] = (x,y)->pie_chart(x,y,MolecularEvolution.sum2one(rand(4)),size = 0.03)
end
tree_draw(newt,draw_labels = false,line_width = 0.5mm, compose_dict = compose_dict)


img = tree_draw(tree)
img |> SVG("imgout.svg",10cm, 10cm)
OR
using Cairo
img |> PDF("imgout.pdf",10cm, 10cm)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/viz/tree_compose.jl#L137-L177" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.tree_polish!-Tuple{Any, Any}' href='#MolecularEvolution.tree_polish!-Tuple{Any, Any}'><span class="jlbinding">MolecularEvolution.tree_polish!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



tree_polish!(newt, models; tol = 10^-4, verbose = 1, topology = true)

Takes a tree and a model function, and optimizes branch lengths and, optionally, topology. Returns final LL. Set `verbose=0` to suppress output. Note: This is not intended for an exhaustive tree search (which requires different heuristics), but rather to polish a tree that is already relatively close to the optimum.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/utils/misc.jl#L160-L165" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.unc2probvec-Tuple{Any}' href='#MolecularEvolution.unc2probvec-Tuple{Any}'><span class="jlbinding">MolecularEvolution.unc2probvec</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
unc2probvec(v)
```


Takes an array of N-1 unbounded values and returns an array of N values that sums to 1. Typically useful for optimizing over categorical probability distributions.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/discrete_models/utils/matrix_helpers.jl#L177-L181" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.univariate_maximize-Tuple{Any, Real, Real, Any, BrentsMethodOpt, Real}' href='#MolecularEvolution.univariate_maximize-Tuple{Any, Real, Real, Any, BrentsMethodOpt, Real}'><span class="jlbinding">MolecularEvolution.univariate_maximize</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
univariate_maximize(f, a::Real, b::Real, transform, optimizer::BrentsMethodOpt, t::Real; ε::Real=sqrt(eps))
```


Maximizes `f(x)` using Brent&#39;s method. See `?brents_method_minimize`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/utils/simple_optim.jl#L198-L202" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.univariate_maximize-Tuple{Any, Real, Real, Any, GoldenSectionOpt, Real}' href='#MolecularEvolution.univariate_maximize-Tuple{Any, Real, Real, Any, GoldenSectionOpt, Real}'><span class="jlbinding">MolecularEvolution.univariate_maximize</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
univariate_maximize(f, a::Real, b::Real, transform, optimizer::GoldenSectionOpt, tol::Real)
```


Maximizes `f(x)` using a Golden Section Search. See `?golden_section_maximize`.

**Examples**

```julia
julia> f(x) = -(x-2)^2
f (generic function with 1 method)

julia> m = univariate_maximize(f, 1, 5, identity, GoldenSectionOpt(), 1e-10)
2.0000000000051843
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/utils/simple_optim.jl#L77-L89" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.univariate_sampler-Tuple{Any, BranchlengthSampler, Any}' href='#MolecularEvolution.univariate_sampler-Tuple{Any, BranchlengthSampler, Any}'><span class="jlbinding">MolecularEvolution.univariate_sampler</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
univariate_sampler(LL, modifier::BranchlengthPeturbation, curr_branchlength)

A MCMC algorithm that draws the next sample of a Markov Chain that approximates the Posterior distrubution over the branchlengths.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/utils/simple_sample.jl#L18-L23" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.values_from_phylo_tree-Tuple{Any, Any}' href='#MolecularEvolution.values_from_phylo_tree-Tuple{Any, Any}'><span class="jlbinding">MolecularEvolution.values_from_phylo_tree</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
values_from_phylo_tree(phylo_tree, key)

Returns a list of values from the given key in the nodes of the phylo_tree, in an order that is somehow compatible with the order the nodes get plotted in.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/viz/phylo_glue.jl#L72-L77" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.weightEM-Tuple{Matrix{Float64}, Any}' href='#MolecularEvolution.weightEM-Tuple{Matrix{Float64}, Any}'><span class="jlbinding">MolecularEvolution.weightEM</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
weightEM(con_lik_matrix::Array{Float64,2}, θ; conc = 0.0, iters = 500)
```


Takes a conditional likelihood matrix (#categories-by-sites) and a starting frequency vector θ (length(θ) = #categories) and optimizes θ (using Expectation Maximization. Maybe.). If conc &gt; 0 then this gives something like variational bayes behavior for LDA. Maybe.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/utils/misc.jl#L40-L45" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.write_fasta-Tuple{String, Vector{String}}' href='#MolecularEvolution.write_fasta-Tuple{String, Vector{String}}'><span class="jlbinding">MolecularEvolution.write_fasta</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
write_fasta(filepath::String, sequences::Vector{String}; seq_names = nothing)
```


Writes a fasta file from a vector of sequences, with optional seq_names.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/utils/fasta_io.jl#L18-L22" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.write_nexus-Tuple{String, FelNode}' href='#MolecularEvolution.write_nexus-Tuple{String, FelNode}'><span class="jlbinding">MolecularEvolution.write_nexus</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
write_nexus(fname::String,tree::FelNode)
```


Writes the tree as a nexus file, suitable for opening in eg. FigTree. Data in the `node_data` dictionary will be converted into annotations. Only tested for simple `node_data` formats and types.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/utils/misc.jl#L281-L287" target="_blank" rel="noreferrer">source</a></Badge>

</details>

