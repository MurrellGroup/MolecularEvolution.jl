
# The MolecularEvolution.jl Framework {#The-MolecularEvolution.jl-Framework}

The organizing principle is that the core algorithms, including Felsenstein&#39;s algorithm, but also a related family of message passing algorithms and inference machinery, are implemented in a way that does not refer to any specific model or even to any particular data type.

## Partitions and BranchModels {#Partitions-and-BranchModels}

A `Partition` is a probabilistic representation of some kind of state. Specifically, it needs to be able to represent P(obs|state) and P(obs,state) when considered as functions of state. So it will typically be able to assign a probability to any possible value of state, and is unnormalized - not required to sum or integrate to 1 over all values of state. As an example, for a discrete state with 4 categories, this could just be a vector of 4 numbers.

For a `Partition` type to be usable by MolecularEvolution.jl, the `combine!` function needs to be implemented. If you have P(obsA|state) and P(obsB|state), then `combine!` calculates P(obsA,obsB|state) under the assumption that obsA and obsB are conditionally independent given state. `MolecularEvolution.jl` tries to avoid allocating memory, so `combine!(dest,src)` places in `dest` the combined `Partition` in `dest`. For a discrete state with 4 categories, this is simply element-wise multiplication of two state vectors.

A `BranchModel` defines how `Partition` distributions evolve along branches. Two functions need to be implemented: `backward!` and `forward!`. We imagine our trees with the root at the top, and `forward!` moves from root to tip, and `backward!` moves from tip to root. `backward!(dest::P,src::P,m::BranchModel,n::FelNode)` takes a src Partition, representing P(obs-below|state-at-bottom-of-branch), and modifies the dest Partition to be P(obs-below|state-at-top-of-branch), where the branch in question is the branch above the `FelNode` n. `forward!` goes in the opposite direction, from P(obs-above,state-at-top-of-branch) to P(obs-above,state-at-bottom-of-branch), with the `Partitions` now, confusingly, representing joint distributions.


![](figures/directions.svg)


## Messages {#Messages}

Nodes on our trees work with messages, where a `message` is a vector of `Partition` structs. This is in case you wish to model multiple different data types on the same tree. Often, all the messages on the tree will just be arrays containing a single `Partition`, but if you&#39;re accessing them you need to remember that they&#39;re in an array!

## Trees {#Trees}

Each node in our tree is a `FelNode` (&quot;Fel&quot; for &quot;Felsenstein&quot;). They point to their parent nodes, and an array of their children, and they store their main vector of `Partition`s, but also cached versions of those from their parents and children, to allow certain message passing schemes. They also have a `branchlength` field, which tells eg. `forward!` and `backward!` how much evolution occurs along the branch above (ie. closer to the root) that node. They also allow for an arbitrary dictionary of `node_data`, in case a model needs any other branch-specific parameters.

The set of algorithms needs to know which model to use for which partition, so the assumption made is that they&#39;ll see an array of models whose order will match the partition array. In general, we might want the models to vary from one branch to another, so the central algorithms take a function that associates a `FelNode->Vector{:<BranchModel}`. In the simpler cases where the model does not vary from branch to branch, or where there is only a single Partition, and thus a single model, the core algorithms have been overloaded to allow you to pass in a single model vector or a single model.

## Algorithms {#Algorithms}

Felsenstein&#39;s algorithm recursively computes, for each node, the probability of all observations below that node, given the state at that node. Felsenstein&#39;s algorithm can be decomposed into the following combination of `backward!` and `combine!` operations:


![](figures/FelsensteinRecursion.svg)


At the root node, we wind up with $P(O_{all}|R)$, where $R$ is the state at the root, and we can compute:

$$P(O_{all}) = \sum_{R} P(O_{all}|R) P(R) \tag{1}$$

::: tip Note

If the root state space would happen to be continuous, then we would instead have:

$$P(O_{all}) = \int_{R} P(O_{all}|R) P(R) \tag{2}$$

:::

After a Felsenstein pass, we can recursively compute, for each node, the probability of all observations above that node, given the state at that node. We call this algorithm `felsenstein_down!`, and it can be decomposed into the following combination of `forward!` and `combine!` operations:


![](figures/FelsensteinDownRecursion.svg)


::: tip Note
- We don&#39;t display how $P(O_{Â}, O_L|C)$ is computed, to avoid cluttering the illustration. It, however, follows from swapping B and C.
  
- When we say &quot;observation above&quot;, we mean the observations in the tree that aren&#39;t below, together with the **root state**.
  

:::

## Technicalities {#Technicalities}

### Scaling constants {#Scaling-constants}

For phylogenetic trees with many nodes or long branch lengths, the likelihood values can become extremely small, leading to numerical underflow in floating-point arithmetic. To maintain numerical stability, we perform computations in the log-domain. This approach is implemented in the following components:
- **[`GaussianPartition`](/api#MolecularEvolution.GaussianPartition) with the field `norm_const`**. Let `part isa GaussianPartition`, then if we integrate the conditional probability of `part` over the real axis, it evaluates to `exp(part.norm_const)`.
  
- **Concrete subtypes of [`DiscretePartition`](/api#MolecularEvolution.DiscretePartition) with the field `scaling`**. Let `part isa DiscretePartition`. The conditional probability, `P`, of state `i` at site `j`, is computed by `P = part.state[i, j] * exp(part.scaling[j])`.
  

### Root state {#Root-state}

In equations (1) &amp; (2) in the above [Algorithms](/framework#Algorithms) section, `parent_message` is the field of a **root**-`FelNode` which represents the quantity `P(R)`. The typical way we specify the root state is to, during initiliazation, pass the desired `parent_message` as the template `Partition` to [`allocate!`](/api#MolecularEvolution.allocate!-Tuple{Any,%20Any}). In addition to the root-`parent_message` being incorporated in computing the likelihood of a tree, it&#39;s also used for downward passes, which includes: [`felsenstein_down!`](/api#MolecularEvolution.felsenstein_down!-Tuple{FelNode,%20Any}), [`branchlength_optim!`](/optimization#MolecularEvolution.branchlength_optim!), [`nni_optim!`](/optimization#MolecularEvolution.nni_optim!), [`marginal_state_dict`](/ancestors#MolecularEvolution.marginal_state_dict).

## Functions {#Functions}
<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.combine!' href='#MolecularEvolution.combine!'><span class="jlbinding">MolecularEvolution.combine!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
combine!(dest::P, src::P) where P<:Partition
```


Combines evidence from two partitions of the same type, storing the result in dest. Note: You should overload this for your own Partititon types.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/e2c752a0b1a5fa2ca25a6c4589b6e4bb030ee730/src/models/discrete_models/discrete_partitions.jl#L185-L190" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.forward!' href='#MolecularEvolution.forward!'><span class="jlbinding">MolecularEvolution.forward!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
forward!(dest::Partition, source::Partition, model::BranchModel, node::FelNode)
```


Propagate the source partition forwards along the branch to the destination partition, under the model. Note: You should overload this for your own BranchModel types.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/e2c752a0b1a5fa2ca25a6c4589b6e4bb030ee730/src/models/discrete_models/pmatrix_models/pmatrix_models.jl#L22-L27" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.backward!' href='#MolecularEvolution.backward!'><span class="jlbinding">MolecularEvolution.backward!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
backward!(dest::Partition, source::Partition, model::BranchModel, node::FelNode)
```


Propagate the source partition backwards along the branch to the destination partition, under the model. Note: You should overload this for your own BranchModel types.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/e2c752a0b1a5fa2ca25a6c4589b6e4bb030ee730/src/models/discrete_models/pmatrix_models/pmatrix_models.jl#L5-L10" target="_blank" rel="noreferrer">source</a></Badge>

</details>

