# # Visualization

# We offer two routes to visualization. The first is using our own plotting routines, built atop Compose.jl. The second converts our trees to Phylo.jl trees, and plots with their Plots.jl recipes. The Compose, Plots, and Phylo dependencies are optional.

# ## Example 1

using MolecularEvolution, Plots, Phylo

#First simulate a tree, and then Brownian motion:
tree = sim_tree(n = 20)
internal_message_init!(tree, GaussianPartition())
bm_model = BrownianMotion(0.0, 0.1)
sample_down!(tree, bm_model)

#We'll add the Gaussian means to the node_data dictionaries
for n in getnodelist(tree)
    n.node_data = Dict(["mu" => n.message[1].mean])
end

#Transducing the mol ev tree to a Phylo.jl tree
phylo_tree = get_phylo_tree(tree)

pl = plot(
    phylo_tree,
    showtips = true,
    tipfont = 6,
    marker_z = "mu",
    markeralpha = 0.5,
    line_z = "mu",
    linecolor = :darkrainbow,
    markersize = 4.0,
    markerstrokewidth = 0,
    margins = 1Plots.cm,
    linewidth = 1.5,
    markercolor = :darkrainbow,
    size = (500, 500),
)

# We also offer `savefig_tweakSVG("simple_plot_example.svg", pl)` for some post-processing tricks that improve the exported trees, like rounding line caps, and `values_from_phylo_tree(phylo_tree,"mu")` which can extract stored quantities in the right order for passing into eg. `markersize` options when plotting.

# For a more comprehensive list of things you can do with Phylo.jl plots, please see [their documentation](https://docs.ecojulia.org/Phylo.jl/stable/man/plotting/).

# ## Drawing trees with `Compose.jl`.

# The `Compose.jl` in-house tree drawing offers extensive flexibility. Here is an example that plots a pie chart representing the marginal probability of each of the 4 possible nucleotides on all nodes on the tree:

using MolecularEvolution, Compose

tree = sim_tree(40, 1000.0, 0.005, mutation_rate = 0.001)
model = DiagonalizedCTMC(reversibleQ(ones(6), ones(4) ./ 4))
internal_message_init!(tree, NucleotidePartition(ones(4) ./ 4, 1))
sample_down!(tree, model)
d = marginal_state_dict(tree, model);
#-
compose_dict = Dict()
for n in getnodelist(tree)
    compose_dict[n] =
        (x, y) -> pie_chart(x, y, d[n][1].state[:, 1], size = 0.02, opacity = 0.75)
end
img = tree_draw(tree,draw_labels = false, line_width = 0.5mm, compose_dict = compose_dict)

# This can then be exported with:

savefig_tweakSVG("piechart_tree.svg",img);

# ## Multiple trees
# Doesn't require `Phylo.jl`. Query trees can be plotted against a reference tree with [`plot_multiple_trees`](@ref). This can be useful, for instance, when we've sampled trees with [`metropolis_sample`](@ref).
using MolecularEvolution, Plots

tree = sim_tree(10, 1, 1)
nodelist = getnodelist(tree); mean = sum([n.branchlength for n in nodelist]) / length(nodelist)
rparams(n::Int) = MolecularEvolution.sum2one(rand(n))
model = DiagonalizedCTMC(reversibleQ(ones(6) ./ (6 * mean), rparams(4)))
internal_message_init!(tree, NucleotidePartition(ones(4) ./ 4, 100))
sample_down!(tree, model)
@time trees, LLs = metropolis_sample(tree, [model], 300, collect_LLs=true);
reference = trees[argmax(LLs)];
# We'll use the maximum a posteriori tree as reference
plot_multiple_trees(trees, reference)
# We can pass in a weight function to fit query trees against `reference` in a weighted least squares fashion with a location and scale parameter.
#=
!!! note
    If we don't want to scale the query trees, we must disable it with `opt_scale = false`.
=#
plot_multiple_trees(
    trees,
    reference,
    y_jitter = 0.05,
    weight_fn = n::FelNode ->
       ifelse(MolecularEvolution.isroot(n) || isleafnode(n), 1.0, 0.0)
)
#jl # We only prefix with `MolecularEvolution` for `isroot`, since we have `Phylo.jl` in our namespace
#=
## Functions

```@docs
get_phylo_tree
values_from_phylo_tree
savefig_tweakSVG
tree_draw
plot_multiple_trees
```
=#
