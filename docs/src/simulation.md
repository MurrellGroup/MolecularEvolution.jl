# Simulation

The two key steps in phylogenetic simulation are 1) simulating the phylogeny itself, and 2) simulating data that evolves over the phylogeny.

## Simulating phylogenies

Author note: while this seems to produce trees with the right "shape", please note that my maths hasn't been extensively checked, and things.

If you just need a simple tree for testing things, then you can just use:

```julia
tree = sim_tree(n=100)
tree_draw(tree, draw_labels = false, canvas_height = 5cm)
```
![](constant_pop_simple.svg)

This has the characteristic "coalescent under constant population size" look.

However, sim_tree is a bit more powerful than this: it aims to simulate branching under a coalescent process with flexible options for how the effective population size, as well as the sampling rate, might change over time. This is important, because the "constant population size" model is quite extreme, and most of the divergence happens in the early internal branches.

A coalescent process runs backwards in time, starting from the most recent tip, and sampling backwards toward the root, coalescing nodes as it goes, and sometimes adding additional sampled tips. With sim_tree, if `nstart` = `add_limit`, then all the tips will be sampled at the same time, and the tree will be ultrametric.

sim_tree has two arguments driving its flexibility. We'll start with "sampling_rate", which controls the rate at which samples are added to the tree. Even under constant effective population size, this can produce interesting behavior.

```julia
for sampling_rate in [5.0, 0.5, 0.05, 0.005]
    tree = sim_tree(100,1000.0,sampling_rate)
    display(tree_draw(tree, draw_labels = false, canvas_height = 5cm))
end
```
![](constant_Ne_samp_rate_5.0.svg)
![](constant_Ne_samp_rate_0.5.svg)
![](constant_Ne_samp_rate_0.05.svg)
![](constant_Ne_samp_rate_0.005.svg)

Above, this rate was just a fixed constant value, but we can also let this be a function. In this example, we'll plot the tree alongside the sampling rate function, as well as the cumulative number of samples through time.

```julia
s(t) = ifelse(0<mod(t/10,10)<1,10.0,0.0)
tree = sim_tree(500,500.0,s)
display(tree_draw(tree, draw_labels = false, canvas_height = 5cm, canvas_width = 23cm))

#Figure out how tall the tree is
root_dists,_ = MolecularEvolution.root2tip_distances(tree)
mrd = maximum(root_dists)
sample_times = mrd .- root_dists
xvals = 0.0:0.1:mrd
display(plot(xvals,s, xflip = true, size = (500,250), xlabel = "time",ylabel = "sampling rate", legend = :none))
display(plot(xvals,x -> sum(x .> sample_times), xflip = true, size = (500,250), xlabel = "time",ylabel = "cumulative samples", legend = :none))
```
![](stepwise_sampling_tree.svg)
![](stepwise_sampling_sampling_rate.svg)
![](stepwise_sampling_cum_samps.svg)

Note how the x axis of these plots is flipped, since the leaf furtherest from the root begins at time=0, and the coalescent runs backwards, from tip to root.

We can also vary the effective population size over time, which adds a different dimension of control. Here is an example showing the shape of a tree under exponential growth:

```julia
n(t) = 100000*exp(-t/10)
tree = sim_tree(100,n,100.0, nstart = 100)
display(tree_draw(tree, draw_labels = false, canvas_height = 7cm, canvas_width = 23cm))

root_dists,_ = MolecularEvolution.root2tip_distances(tree)
plot(0.0:0.1:maximum(root_dists),n, xflip = true, size = (500,250), xlabel = "time",ylabel = "effective population size", legend = :none)
```
![](exp_growth_tree.svg)
![](exp_growth_popsize.svg)

Logistic growth, with a relatively low sampling rate, provides a reasonable model of an emerging virus that was only sampled later in its growth trajectory, such as HIV.

```julia
n(t) = 10000/(1+exp(t-10))
tree = sim_tree(100,n,20.0)
display(tree_draw(tree, draw_labels = false, canvas_height = 7cm, canvas_width = 23cm))

root_dists,_ = MolecularEvolution.root2tip_distances(tree)
display(plot(0.0:0.1:maximum(root_dists),n, xflip = true, size = (500,250), xlabel = "time",ylabel = "effective population size", legend = :none))

mrd = maximum(root_dists)
sample_times = mrd .- root_dists
plot(0.0:0.1:mrd,x -> sum(x .> sample_times), xflip = true, size = (500,250), xlabel = "time",ylabel = "cumulative samples", legend = :none)
```
![](logistic_growth_tree.svg)
![](logistic_growth_popsize.svg)
![](logistic_growth_cum_samps.svg)

How about a virus with a seasonally varying effective population size, where sampling is proportional to case counts? Between seasons, the effective population size gets so low that the next seasons clade arises from a single lineage in the previous season.

```julia
n(t) = exp(sin(t/10) * 2.0 + 4)
s(t) = n(t)/100
tree = sim_tree(500,n,s)
display(tree_draw(tree, draw_labels = false))


root_dists,_ = MolecularEvolution.root2tip_distances(tree)
display(plot(0.0:0.1:maximum(root_dists),n, xflip = true, size = (500,250), xlabel = "time",ylabel = "effective population size", legend = :none))

mrd = maximum(root_dists)
sample_times = mrd .- root_dists
plot(0.0:0.1:mrd,x -> sum(x .> sample_times), xflip = true, size = (500,250), xlabel = "time",ylabel = "cumulative samples", legend = :none)
```
![](seasonal_tree.svg)
![](seasonal_popsize.svg)
![](seasonal_cum_samps.svg)

Finally, the `mutation_rate` argument multiplicatively scales the branch lengths.

## Simulating evolution over phylogenies

Coming soon.