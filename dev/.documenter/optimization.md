
# Optimization {#Optimization}

There are two distinct kinds of optimization: &quot;global&quot; model parameters, and then tree branchlengths and topology. These are kept distinct because we can use algorithmic tricks to dramatically improve the performance of the latter.

The example below will set up and optimize a [&quot;Generalized Time Reversible&quot; nucleotide substitution model](https://en.wikipedia.org/wiki/Substitution_model), where there are 6 rate parameters that govern the symmetric part of a rate matrix, and 4 nucleotide frequencies (that sum to 1, so only 3 underlying parameters).

::: tip Note

For the Bayesian counterpart of this page, we refer you to [Set up a Bayesian model sampler](/generated/update#Set-up-a-Bayesian-model-sampler) and  [Multiple trees](/generated/viz#Multiple-trees).

:::

## Optimizing model parameters {#Optimizing-model-parameters}

We first need to construct an objective function. A very common use case involves parameterizing a rate matrix (along with all the constraints this entails) from a flat parameter vector. `reversibleQ` can be convenient here, which takes a vector of parameters and equilibrium frequencies and returns a reversible rate matrix. The parameters are the upper triangle (excluding the diagonal) of the rate matrix:

```julia
reversibleQ(1:6,ones(4))
```


```
4×4 Matrix{Float64}:
 -6.0    1.0    2.0    3.0
  1.0  -10.0    4.0    5.0
  2.0    4.0  -12.0    6.0
  3.0    5.0    6.0  -14.0
```


...and the equilibrium frequencies are multiplied column-wise:

```julia
reversibleQ(ones(6),[0.1,0.2,0.3,0.4])
```


```
4×4 Matrix{Float64}:
 -0.9   0.2   0.3   0.4
  0.1  -0.8   0.3   0.4
  0.1   0.2  -0.7   0.4
  0.1   0.2   0.3  -0.6
```


Another convenient trick is to be able to parameterize a vector of positive frequencies that sum to 1, using N-1 unconstrained parameters. `unc2probvec` can help:

```julia
unc2probvec(zeros(3))
```


```
4-element Vector{Float64}:
 0.25
 0.25
 0.25
 0.25
```


[`ParameterHandling.jl`](https://github.com/invenia/ParameterHandling.jl) provides a convenient framework for managing collections of parameters in a way that plays with much of the Julia optimization ecosystem, and we recommend its use. Here we&#39;ll use `ParameterHandling` and [`NLopt`](https://github.com/JuliaOpt/NLopt.jl).

First, we&#39;ll load in some example nucleotide data:

```julia
using MolecularEvolution, FASTX, ParameterHandling, NLopt

#Read in seqs and tree, and populate the three  NucleotidePartitions
seqnames, seqs = read_fasta("Data/MusNuc_IGHV.fasta")
tree = read_newick_tree("Data/MusNuc_IGHV.tre")
initial_partition = NucleotidePartition(length(seqs[1]))
populate_tree!(tree,initial_partition,seqnames,seqs)
```


Then we set up the model parameters, and the objective function:

```julia
#Named tuple of parameters, with initial values and constraints (from ParameterHandling.jl)
initial_params = (
        rates=positive(ones(6)), #rates must be non-negative
        pi=zeros(3) #will be transformed into 4 eq freqs
)
flat_initial_params, unflatten = value_flatten(initial_params) #See ParameterHandling.jl docs
num_params = length(flat_initial_params)

#Set up a function that builds a model from these parameters
function build_model_vec(params)
    pi = unc2probvec(params.pi)
    return DiagonalizedCTMC(reversibleQ(params.rates,pi))
end

#Set up the function to be *minimized*
function objective(params::NamedTuple; tree = tree)
    #In this example, we are optimizing the nuc equilibrium freqs
    #We'll also assume that the starting frequencies (at the root of the tree) are the eq freqs
    tree.parent_message[1].state .= unc2probvec(params.pi)
    return -log_likelihood!(tree,build_model_vec(params)) #Note, negative of LL, because minimization
end
```


Then we&#39;ll set up an optimizer from `NLOpt`. See [this discussion](https://discourse.julialang.org/t/optim-what-optimiser-is-best-if-your-gradient-computation-is-slow/5487/12) and [this exploration](https://github.com/SciML/DiffEqParamEstim.jl/blob/master/test/lorenz_true_test.jl) of optimizers.

```julia
opt = Opt(:LN_BOBYQA, num_params)
#Note: NLopt requires a function that returns a gradient, even for gradient free methods, hence (x,y)->...
min_objective!(opt, (x,y) -> (objective ∘ unflatten)(x)) #See ParameterHandling.jl docs for objective ∘ unflatten explanation
#Some bounds (which will be in the transformed domain) to prevent searching numerically silly bits of parameter space:
lower_bounds!(opt, [-10.0 for i in 1:num_params])
upper_bounds!(opt, [10.0 for i in 1:num_params])
xtol_rel!(opt, 1e-12)
_,mini,_ = NLopt.optimize(opt, flat_initial_params)
final_params = unflatten(mini)

optimized_model = build_model_vec(final_params)
println("Opt LL:",log_likelihood!(tree,optimized_model))
```


```
Opt LL:-3783.226756522292
```


We can view the optimized parameter values:

```
println("Rates: ", round.(final_params.rates,sigdigits = 4))
println("Pi:", round.(unc2probvec(final_params.pi),sigdigits = 4))
```


```
Rates: [1.124, 2.102, 1.075, 0.9802, 1.605, 0.5536]
Pi:[0.2796, 0.2192, 0.235, 0.2662]
```


Or the entire optimized rate matrix:

```julia
matrix_for_display(optimized_model.Q,['A','C','G','T'])
```


```
Opt LL:-3783.226756522292
5×5 Matrix{Any}:
 ""     'A'        'C'        'G'        'T'
 'A'  -1.02672    0.246386   0.494024   0.286309
 'C'   0.314289  -0.971998   0.23034    0.427368
 'G'   0.587774   0.214842  -0.950007   0.147391
 'T'   0.300663   0.35183    0.130093  -0.782586
```


## Optimizing the tree topology and branch lengths {#Optimizing-the-tree-topology-and-branch-lengths}

With a tree and a model, we can also optimize the branch lengths and search, by [nearest neighbour interchange](https://en.wikipedia.org/wiki/Tree_rearrangement) for changes to the tree that improve the likelihood. Individually, these are performed by `nni_optim!` and `branchlength_optim!`, which need to have `felsenstein!` and `felsenstein_down!` called beforehand, but this is all bundled into:

```julia
tree_polish!(tree, optimized_model)
```


```
LL: -3783.226756522292
LL: -3782.345818028071
LL: -3782.3231632207567
LL: -3782.3211724011044
LL: -3782.321068684831
LL: -3782.3210622627776
```


And just to convince you this works, we can perturb the branch lengths, and see how the likelihood improves:

```julia
for n in getnodelist(tree)
    n.branchlength *= (rand()+0.5)
end
tree_polish!(tree, optimzed_model)
```


```
LL: -3805.4140940138795
LL: -3782.884883999107
LL: -3782.351780962518
LL: -3782.322906364547
LL: -3782.321183009534
LL: -3782.3210398963506
LL: -3782.3210271696703
```


::: warning Warning

`tree_polish!` probably won&#39;t find a good tree from a completely start. Different tree search heuristics are required for that.

:::

## Functions {#Functions}
<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.reversibleQ' href='#MolecularEvolution.reversibleQ'><span class="jlbinding">MolecularEvolution.reversibleQ</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
reversibleQ(param_vec,eq_freqs)
```


Takes a vector of parameters and equilibrium frequencies and returns a reversible rate matrix. The parameters are the upper triangle of the rate matrix, with the diagonal elements omitted, and the equilibrium frequencies are multiplied column-wise.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/discrete_models/utils/matrix_helpers.jl#L122-L128" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.unc2probvec' href='#MolecularEvolution.unc2probvec'><span class="jlbinding">MolecularEvolution.unc2probvec</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
unc2probvec(v)
```


Takes an array of N-1 unbounded values and returns an array of N values that sums to 1. Typically useful for optimizing over categorical probability distributions.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/models/discrete_models/utils/matrix_helpers.jl#L177-L181" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.branchlength_optim!' href='#MolecularEvolution.branchlength_optim!'><span class="jlbinding">MolecularEvolution.branchlength_optim!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



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
<summary><a id='MolecularEvolution.nni_optim!' href='#MolecularEvolution.nni_optim!'><span class="jlbinding">MolecularEvolution.nni_optim!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



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
<summary><a id='MolecularEvolution.root_optim!' href='#MolecularEvolution.root_optim!'><span class="jlbinding">MolecularEvolution.root_optim!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



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
<summary><a id='MolecularEvolution.tree_polish!' href='#MolecularEvolution.tree_polish!'><span class="jlbinding">MolecularEvolution.tree_polish!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



tree_polish!(newt, models; tol = 10^-4, verbose = 1, topology = true)

Takes a tree and a model function, and optimizes branch lengths and, optionally, topology. Returns final LL. Set `verbose=0` to suppress output. Note: This is not intended for an exhaustive tree search (which requires different heuristics), but rather to polish a tree that is already relatively close to the optimum.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/utils/misc.jl#L160-L165" target="_blank" rel="noreferrer">source</a></Badge>

</details>

