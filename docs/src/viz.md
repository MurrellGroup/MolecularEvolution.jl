# Visualization

We offer two routes to visualization. The first is using our own plotting routines, built atop Compose.jl. The second converts our trees to Phylo.jl trees, and plots with their Plots.jl recipes. The Compose, Plots, and Phylo dependencies are optional.

## Example 1

```julia
using MolecularEvolution, Plots, Phylo

#First simulate a tree, and then Brownian motion:
tree = sim_tree(n=20)
internal_message_init!(tree, GaussianPartition())
bm_model = BrownianMotion(0.0,0.1)
sample_down!(tree, bm_model)

#We'll add the Gaussian means to the node_data dictionaries
for n in getnodelist(tree)
    n.node_data = Dict(["mu"=>n.message[1].mean])
end

#Transducing the mol ev tree to a Phylo.jl tree
phylo_tree = get_phylo_tree(tree)

pl = plot(phylo_tree,
    showtips = true, tipfont = 6, marker_z = "mu", markeralpha = 0.5, line_z = "mu", linecolor = :darkrainbow, 
    markersize = 4.0, markerstrokewidth = 0,margins = 1Plots.cm,
    linewidth = 1.5, markercolor = :darkrainbow, size = (500, 500))
```

![](figures/simple_plot_example.svg)

We also offer `savefig_tweakSVG("simple_plot_example.svg", plot = pl)` for some post-processing tricks that improve the exported trees, like rounding line caps, and `values_from_phylo_tree(phylo_tree,"mu")` which can extract stored quantities in the right order for passing into eg. `markersize` options when plotting.

For a more comprehensive list of things you can do with Phylo.jl plots, please see [their documentation](https://docs.ecojulia.org/Phylo.jl/stable/man/plotting/).

## Coming soon.
Examples using Compose dependent plots.