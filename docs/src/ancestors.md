# Ancestral Reconstruction

Given a phylogeny, and observations on some set of leaf nodes, "ancestral reconstruction" describes a family of approaches for inferring the state of the ancestors, or the distribution over possible states of ancestors.

## Examples

```julia
using MolecularEvolution

#Simulate a small tree, with Brownian motion over it
tree = sim_tree(n=10)
internal_message_init!(tree, GaussianPartition())
bm_model = BrownianMotion(0.0,0.1)
sample_down!(tree, bm_model)

r(x) = round(x,sigdigits = 3)
println("Leaf values:")
for n in getleaflist(tree)
    println(n.name," : ",r(n.message[1].mean))
end

d = marginal_state_dict(tree,bm_model)
println("Inferred internal means (±95% intervals):")
for n in getnonleaflist(tree)
    m,s = d[n][1].mean,sqrt(d[n][1].var)
    println(r(m), "±", r(1.96*s), " - true value: ",r(n.message[1].mean))
end
```
```
Leaf values:
tax8 : -1.03
tax1 : -1.15
tax9 : -1.67
tax10 : -0.112
tax6 : -0.0183
tax2 : -0.0574
tax3 : 0.207
tax5 : 0.0021
tax4 : 0.634
tax7 : 0.544
Inferred internal means (±95% intervals):
-0.485±0.815 - true value: -0.587
-1.17±0.556 - true value: -1.37
-1.1±0.256 - true value: -1.09
0.116±0.45 - true value: 0.21
0.0275±0.35 - true value: -0.035
0.0216±0.283 - true value: 0.0177
0.0459±0.13 - true value: 0.0485
0.0532±0.122 - true value: 0.075
0.571±0.147 - true value: 0.589
```

We can also find the values of the state for each node that are jointly most likely (in the case of Brownian motion, these just happen to be the same as the marginal means, but that isn't the case for other models):

```julia
d = max_joint_state_dict(tree,bm_model)
println("Inferred most likely (jointly) internal values:")
for n in getnonleaflist(tree)
    m = d[n][1].mean
    println(r(m), " - true value: ",r(n.message[1].mean))
end
```
```
Inferred most likely (jointly) internal values:
-0.485 - true value: -0.587
-1.17 - true value: -1.37
-1.1 - true value: -1.09
0.116 - true value: 0.21
0.0275 - true value: -0.035
0.0216 - true value: 0.0177
0.0459 - true value: 0.0485
0.0532 - true value: 0.075
0.571 - true value: 0.589
```

And we can sample internal states under our model, but conditioned on the leaf observations:

```julia
d = endpoint_conditioned_sample_state_dict(tree,bm_model)
println("Sampled states, conditioned on observed leaves:")
for n in getnonleaflist(tree)
    m = d[n][1].mean
    println(r(m), " - true value: ",r(n.message[1].mean))
end
```
```
Sampled states, conditioned on observed leaves:
-0.784 - true value: -0.587
-1.3 - true value: -1.37
-1.13 - true value: -1.09
-0.155 - true value: 0.21
0.0118 - true value: -0.035
0.0305 - true value: 0.0177
0.0913 - true value: 0.0485
0.0542 - true value: 0.075
0.498 - true value: 0.589
```

## Functions

```@docs
marginal_state_dict
max_joint_state_dict
endpoint_conditioned_sample_state_dict
```

