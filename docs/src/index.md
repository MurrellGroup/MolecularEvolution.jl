```@meta
CurrentModule = MolecularEvolution
```

# MolecularEvolution

Documentation for [MolecularEvolution](https://github.com/MurrellGroup/MolecularEvolution.jl).

### A Julia package for the flexible development of phylogenetic models.

MolecularEvolution.jl exploits Julia's multiple dispatch, implementing a fully generic suite of likelihood calculations, branchlength optimization, topology optimization, and ancestral inference. Users can construct trees using already-defined data types and models. But users can define probability distributions over their own data types, and specify the behavior of these under their own model types, and can mix and match different models on the same phylogeny.

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

```@index
```

```@autodocs
Modules = [MolecularEvolution]
```
