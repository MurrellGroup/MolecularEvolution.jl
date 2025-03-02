# Intro

MolecularEvolution.jl exploits Julia's multiple dispatch, implementing a fully generic suite of likelihood calculations, branchlength optimization, topology optimization, and ancestral inference. Users can construct trees using already-defined data types and models. But users can define probability distributions over their own data types, and specify the behavior of these under their own model types, and can mix and match different models on the same phylogeny.

If the behavior you need is not already available in `MolecularEvolution.jl`:
- If you have a new data type:
  - A `Partition` type that represents the uncertainty over your state. 
  - `combine!()` that merges evidence from two `Partition`s.
- If you have a new model:
  - A `BranchModel` type that stores your model parameters.
  - `forward!()` that evolves state distributions over branches, in the root-to-tip direction.
  - `backward!()` that reverse-evolves state distributions over branches, in the tip-to-root direction.

And then sampling, likelihood calculations, branch-length optimization, ancestral reconstruction, etc should be available for your new data or model.

### Design principles
In order of importance, we aim for the following:
- Flexibility and generality
  - Where possible, we avoid design decisions that limit the development of new models, or make it harder to develop new models.
  - We do not sacrifice flexibility for performance.
- Scalability
  - Analyses implemented using `MolecularEvolution.jl` should scale to large, real-world datasets.
- Performance
  - While the above take precedence over speed, it should be possible to optimize your `Partition`, `combine!()`, `BranchModel`, `forward!()` and `backward!()` functions to obtain competative runtimes.

### Authors:
Venkatesh Kumar and Ben Murrell, with additional contributions by Sanjay Mohan, Alec Pankow, Hassan Sadiq, and Kenta Sato.

### Quick example: Likelihood calculations under phylogenetic Brownian motion:

```julia
using MolecularEvolution, Plots

#First simulate a tree, using a coalescent process
tree = sim_tree(n=200)
internal_message_init!(tree, GaussianPartition())
#Simulate brownian motion over the tree
bm_model = BrownianMotion(0.0,1.0)
sample_down!(tree, bm_model)
#And plot the log likelihood as a function of the parameter value
ll(x) = log_likelihood!(tree,BrownianMotion(0.0,x))
plot(0.7:0.001:1.6,ll, xlabel = "variance per unit time", ylabel = "log likelihood")
```
![](figures/quick_example.svg)
