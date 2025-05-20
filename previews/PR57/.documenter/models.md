
# Models {#Models}

We offer a catalogue of frequently used models that are already integrated with the framework and ready to be used. We maintain that if a model you&#39;d like to use is not included in the list, you can swiftly define one yourself and leverage our framework nonetheless.

## Discrete state models {#Discrete-state-models}

Here we account for a typical set-up for a discrete state `Partition`:
<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.DiscretePartition-models' href='#MolecularEvolution.DiscretePartition-models'><span class="jlbinding">MolecularEvolution.DiscretePartition</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
abstract type DiscretePartition <: MultiSitePartition end
```


Represents a `states`-by-`sites` matrix of probabilities. The following fields are loosely required:
- `state`: A matrix of probabilities that are site-wise normalized.
  
- `states`: The number of states.
  
- `sites`: The number of sites.
  
- `scaling`: A vector of log-domain probability scaling factors, one per site.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/e2c752a0b1a5fa2ca25a6c4589b6e4bb030ee730/src/models/discrete_models/discrete_partitions.jl#L5-L13" target="_blank" rel="noreferrer">source</a></Badge>

</details>


And here&#39;s a list of simple concrete subtypes of `DiscretePartition`:
- [`CodonPartition`](/api#MolecularEvolution.CodonPartition)
  
- [`CustomDiscretePartition`](/api#MolecularEvolution.CustomDiscretePartition)
  
- [`NucleotidePartition`](/api#MolecularEvolution.NucleotidePartition)
  
- [`GappyNucleotidePartition`](/api#MolecularEvolution.GappyNucleotidePartition)
  
- [`AminoAcidPartition`](/api#MolecularEvolution.AminoAcidPartition)
  
- [`GappyAminoAcidPartition`](/api#MolecularEvolution.GappyAminoAcidPartition)
  

And then there are two typical `BranchModel`s that will cooperate with this `Partition`:
- [`GeneralCTMC`](/api#MolecularEvolution.GeneralCTMC)
  
- [`DiagonalizedCTMC`](/api#MolecularEvolution.DiagonalizedCTMC)
  

::: tip Note

The two above can be regarded as special cases of the more general [`PModel`](/api#MolecularEvolution.PModel), which just represents a P matrix.

:::

A typical way of constructing your Q matrix in our ecosystem is by:
<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.reversibleQ-models' href='#MolecularEvolution.reversibleQ-models'><span class="jlbinding">MolecularEvolution.reversibleQ</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
reversibleQ(param_vec,eq_freqs)
```


Takes a vector of parameters and equilibrium frequencies and returns a reversible rate matrix. The parameters are the upper triangle of the rate matrix, with the diagonal elements omitted, and the equilibrium frequencies are multiplied column-wise.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/e2c752a0b1a5fa2ca25a6c4589b6e4bb030ee730/src/models/discrete_models/utils/matrix_helpers.jl#L122-L128" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.nonreversibleQ-models' href='#MolecularEvolution.nonreversibleQ-models'><span class="jlbinding">MolecularEvolution.nonreversibleQ</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
nonreversibleQ(param_vec)
```


Takes a vector of parameters and returns a nonreversible rate matrix.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/e2c752a0b1a5fa2ca25a6c4589b6e4bb030ee730/src/models/discrete_models/utils/matrix_helpers.jl#L151-L155" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### Codon models {#Codon-models}
<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.CodonPartition-models' href='#MolecularEvolution.CodonPartition-models'><span class="jlbinding">MolecularEvolution.CodonPartition</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



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


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/e2c752a0b1a5fa2ca25a6c4589b6e4bb030ee730/src/models/discrete_models/codon_models.jl#L349-L357" target="_blank" rel="noreferrer">source</a></Badge>

</details>


We offer constructors for the following Q matrix parameterizations:
- `MolecularEvolution.MG94_F3x4` - for example usage, see [Example 3: FUBAR](/examples#Example-3:-FUBAR)
  
- `MolecularEvolution.MG94_F61`
  
- `MolecularEvolution.HB98_F61`
  
- `MolecularEvolution.HB98_AAfit`
  

Use help mode, `?`, in the REPL to find out more about the above.

### Miscellaneous models {#Miscellaneous-models}
<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.InterpolatedDiscreteModel-models' href='#MolecularEvolution.InterpolatedDiscreteModel-models'><span class="jlbinding">MolecularEvolution.InterpolatedDiscreteModel</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



**Constructors**

```julia
InterpolatedDiscreteModel(siz::Int64, generator, tvec::Vector{Float64})
InterpolatedDiscreteModel(Pvec::Array{Float64,3}, tvec::Vector{Float64})
```


`generator` is a function that takes a time value `t` and returns a P matrix.

**Description**

Stores a number (`siz`) of P matrices, and the time values to which they correspond. For a requested t, the returned P matrix is (element-wise linearly) interpolated between it&#39;s two neighbours.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/e2c752a0b1a5fa2ca25a6c4589b6e4bb030ee730/src/models/discrete_models/interpolated_discrete_model.jl#L5-L15" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.PiQ-models' href='#MolecularEvolution.PiQ-models'><span class="jlbinding">MolecularEvolution.PiQ</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



**Constructors**

```julia
PiQ(r::Float64,pi::Vector{Float64}; normalize=false)
PiQ(pi::Vector{Float64}; normalize=false)
```


**Description**

The F81 substitution model, but for general dimensions. https://www.diva-portal.org/smash/get/diva2:1878793/FULLTEXT01.pdf


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/e2c752a0b1a5fa2ca25a6c4589b6e4bb030ee730/src/models/discrete_models/PiQ.jl#L4-L13" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Continuous models {#Continuous-models}

The partition of interest is:
<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.GaussianPartition-models' href='#MolecularEvolution.GaussianPartition-models'><span class="jlbinding">MolecularEvolution.GaussianPartition</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



**Constructors**

```julia
GaussianPartition(mean::Float64, var::Float64, norm_const::Float64)
GaussianPartition(mean::Float64, var::Float64) # norm_const defaults to 0.0
GaussianPartition() # mean, var, norm_const default to 0.0, 1.0, 0.0 respectively
```


**Description**

A partition representing a (not necessarily normalized) Gaussian distribution. `norm_const` is the log-domain normalization constant.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/e2c752a0b1a5fa2ca25a6c4589b6e4bb030ee730/src/models/continuous_models/gaussian_partition.jl#L3-L12" target="_blank" rel="noreferrer">source</a></Badge>

</details>


And then there are two `BranchModel`s which are compatible with the above partition, namely:
<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.BrownianMotion-models' href='#MolecularEvolution.BrownianMotion-models'><span class="jlbinding">MolecularEvolution.BrownianMotion</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
BrownianMotion(mean_drift::Float64, var_drift::Float64)
```


A 1D continuous Brownian motion model with mean drift and variance drift.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/e2c752a0b1a5fa2ca25a6c4589b6e4bb030ee730/src/models/continuous_models/brownian_motion.jl#L2-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.ZeroDriftOrnsteinUhlenbeck-models' href='#MolecularEvolution.ZeroDriftOrnsteinUhlenbeck-models'><span class="jlbinding">MolecularEvolution.ZeroDriftOrnsteinUhlenbeck</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
function ZeroDriftOrnsteinUhlenbeck(
    mean::Float64,
    volatility::Float64,
    reversion::Float64,
)
```


A 1D continuous Ornstein-Uhlenbeck process with mean, volatility, and reversion.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/e2c752a0b1a5fa2ca25a6c4589b6e4bb030ee730/src/models/continuous_models/ornstein_uhlenbeck.jl#L3-L11" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Compound models {#Compound-models}

### Branch-wise mixture {#Branch-wise-mixture}
<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.BWMModel-models' href='#MolecularEvolution.BWMModel-models'><span class="jlbinding">MolecularEvolution.BWMModel</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



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


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/e2c752a0b1a5fa2ca25a6c4589b6e4bb030ee730/src/models/compound_models/bwm.jl#L2-L13" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### CAT {#CAT}
<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.CATModel-models' href='#MolecularEvolution.CATModel-models'><span class="jlbinding">MolecularEvolution.CATModel</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
CATModel(models::Vector{<:BranchModel})
```


CAT is something where you split the sites up, and assign each site to a different model (whose &quot;data&quot; gets stored in a contiguous block of memory).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/e2c752a0b1a5fa2ca25a6c4589b6e4bb030ee730/src/models/compound_models/cat.jl#L2-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.CATPartition-models' href='#MolecularEvolution.CATPartition-models'><span class="jlbinding">MolecularEvolution.CATPartition</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



**Constructors**

```julia
CATPartition(part_inds::Vector{Vector{Int}})
CATPartition(part_inds::Vector{Vector{Int}}, parts::Vector{PType})
```


**Description**

A partition for the [`CATModel`](/api#MolecularEvolution.CATModel).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/e2c752a0b1a5fa2ca25a6c4589b6e4bb030ee730/src/models/compound_models/cat.jl#L15-L23" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### Covarion {#Covarion}
<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.CovarionModel-models' href='#MolecularEvolution.CovarionModel-models'><span class="jlbinding">MolecularEvolution.CovarionModel</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



**Constructors**

```julia
CovarionModel(models::Vector{<:DiscreteStateModel}, inter_transition_rates::Matrix{Float64})
CovarionModel(models::Vector{<:DiscreteStateModel}, inter_transition_rate::Float64)
```


**Description**

The covarion model.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/e2c752a0b1a5fa2ca25a6c4589b6e4bb030ee730/src/models/compound_models/covarion.jl#L2-L10" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.CovarionPartition-models' href='#MolecularEvolution.CovarionPartition-models'><span class="jlbinding">MolecularEvolution.CovarionPartition</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
CovarionPartition(states,sites,models,t)
```


A partition for the [`CovarionModel`](/api#MolecularEvolution.CovarionModel).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/e2c752a0b1a5fa2ca25a6c4589b6e4bb030ee730/src/models/compound_models/covarion.jl#L47-L51" target="_blank" rel="noreferrer">source</a></Badge>

</details>


### Site-wise mixture {#Site-wise-mixture}

See [Example 2: GTR+Gamma](/examples#Example-2:-GTRGamma) for example usage.
<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.SWMModel-models' href='#MolecularEvolution.SWMModel-models'><span class="jlbinding">MolecularEvolution.SWMModel</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



**Constructors**

```julia
SWMModel(models::Vector{<:BranchModel})
SWMModel(model::M, rs::Vector{Float64}) where {M <: BranchModel}
```


**Description**

A site-wise mixture model, for site-to-site &quot;random effects&quot; rate variation.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/e2c752a0b1a5fa2ca25a6c4589b6e4bb030ee730/src/models/compound_models/swm.jl#L17-L25" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.SWMPartition-models' href='#MolecularEvolution.SWMPartition-models'><span class="jlbinding">MolecularEvolution.SWMPartition</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



**Constructors**

```julia
SWMPartition(parts::Vector{PType}) where {PType <: MultiSitePartition}
SWMPartition(part::PType, n_parts::Int) where {PType <: MultiSitePartition}
SWMPartition(parts::Vector{PType}, weights::Vector{Float64}, sites::Int, states::Int, models::Int) where {PType <: MultiSitePartition}
```


**Description**

A site-wise mixture partition for the [`SWMModel`](/api#MolecularEvolution.SWMModel).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/e2c752a0b1a5fa2ca25a6c4589b6e4bb030ee730/src/models/compound_models/swm.jl#L46-L55" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Lazy models {#Lazy-models}

### LazyPartition {#LazyPartition}
<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.LazyPartition-models' href='#MolecularEvolution.LazyPartition-models'><span class="jlbinding">MolecularEvolution.LazyPartition</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



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
    
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/e2c752a0b1a5fa2ca25a6c4589b6e4bb030ee730/src/models/lazy_models/lazy_partition.jl#L2-L25" target="_blank" rel="noreferrer">source</a></Badge>

</details>


#### Examples {#Examples}

##### Example 1: Initializing for an upward pass {#Example-1:-Initializing-for-an-upward-pass}

Now, we show how to wrap the `CodonPartition`s from [Example 3: FUBAR](/examples#Example-3:-FUBAR) with `LazyPartition`:

You simply go from initializing messages like this:

```julia
initial_partition = CodonPartition(Int64(length(seqs[1])/3))
initial_partition.state .= eq_freqs
populate_tree!(tree,initial_partition,seqnames,seqs)
```


To this

```julia
initial_partition = CodonPartition(Int64(length(seqs[1])/3))
initial_partition.state .= eq_freqs
lazy_initial_partition = LazyPartition{CodonPartition}()
populate_tree!(tree,lazy_initial_partition,seqnames,seqs)
lazyprep!(tree, initial_partition)
```


By this slight modification, we go from initializing and using 554 partitions to 6 during the subsequent `log_likelihood!` and `felsenstein!` calls. There is no significant decrease in performance recorded from this switch.

##### Example 2: Initializing for a downward pass {#Example-2:-Initializing-for-a-downward-pass}

Now, we show how to wrap the `GaussianPartition`s from [Quick example: Likelihood calculations under phylogenetic Brownian motion:](/intro#Quick-example:-Likelihood-calculations-under-phylogenetic-Brownian-motion:) with `LazyPartition`:

You simply go from initializing messages like this:

```julia
internal_message_init!(tree, GaussianPartition())
```


To this (technically we only add 1 LOC)

```julia
initial_partition = GaussianPartition()
lazy_initial_partition = LazyPartition{GaussianPartition}()
internal_message_init!(tree, lazy_initial_partition)
lazyprep!(tree, initial_partition, direction=LazyDown(isleafnode))
```


::: tip Note

Now, we provided a direction for `lazyprep!`. The direction is an instance of `LazyDown`, which was initialized with the `isleafnode` function. The function `isleafnode` dictates if a node saves its sampled observation after a down pass. If you use `direction=LazyDown()`, every node saves its observation.

:::

#### Surrounding LazyPartition {#Surrounding-LazyPartition}
<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.lazyprep!-models' href='#MolecularEvolution.lazyprep!-models'><span class="jlbinding">MolecularEvolution.lazyprep!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
lazyprep!(tree::FelNode, initial_message::Vector{<:Partition}; partition_list = 1:length(tree.message), direction::LazyDirection = LazyUp())
```


Extra, intermediate step of tree preparations between initializing messages across the tree and calling message passing algorithms with `LazyPartition`.
1. Perform a `lazysort!` on tree to obtain the optimal tree for a lazy `felsenstein!` prop, or a `sample_down!`.
  
2. Fix `tree.parent_message` to an initial message.
  
3. Preallocate sufficiently many inner partitions needed for a `felsenstein!` prop, or a `sample_down!`.
  
4. Specialized preparations based on the direction of the operations (`forward!`, `backward!`). `LazyDown` or `LazyUp`.
  

See also `LazyDown`, `LazyUp`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/e2c752a0b1a5fa2ca25a6c4589b6e4bb030ee730/src/models/lazy_models/lazy_partition.jl#L200-L210" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.LazyUp-models' href='#MolecularEvolution.LazyUp-models'><span class="jlbinding">MolecularEvolution.LazyUp</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



**Constructor**

```
LazyUp()
```


**Description**

Indicate that we want to do an upward pass, e.g. `felsenstein!`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/e2c752a0b1a5fa2ca25a6c4589b6e4bb030ee730/src/models/lazy_models/lazy_partition.jl#L42-L48" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.LazyDown-models' href='#MolecularEvolution.LazyDown-models'><span class="jlbinding">MolecularEvolution.LazyDown</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



**Constructors**

```
LazyDown(stores_obs)
LazyDown() = LazyDown(x::FelNode -> true)
```


**Description**

Indicate that we want to do a downward pass, e.g. `sample_down!`. The function passed to the constructor takes a `node::FelNode` as input and returns a `Bool` that decides if `node` stores its observations.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/e2c752a0b1a5fa2ca25a6c4589b6e4bb030ee730/src/models/lazy_models/lazy_partition.jl#L51-L59" target="_blank" rel="noreferrer">source</a></Badge>

</details>

