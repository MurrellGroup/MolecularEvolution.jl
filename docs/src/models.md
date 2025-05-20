# Models

We offer a catalogue of frequently used models that are already integrated with the framework and ready to be used.
We maintain that if a model you'd like to use is not included in the list, you can swiftly define one yourself and
leverage our framework nonetheless.

## Discrete state models

Here we account for a typical set-up for a discrete state `Partition`:
```@docs; canonical=false
DiscretePartition
```

And here's a list of simple concrete subtypes of `DiscretePartition`:
- [`CodonPartition`](@ref)
- [`CustomDiscretePartition`](@ref)
- [`NucleotidePartition`](@ref)
- [`GappyNucleotidePartition`](@ref)
- [`AminoAcidPartition`](@ref)
- [`GappyAminoAcidPartition`](@ref)

And then there are two typical `BranchModel`s that will cooperate with this `Partition`:
- [`GeneralCTMC`](@ref)
- [`DiagonalizedCTMC`](@ref)
!!! note
    The two above can be regarded as special cases of the more general [`PModel`](@ref), which just represents a P matrix.

A typical way of constructing your Q matrix in our ecosystem is by:
```@docs; canonical=false
reversibleQ
nonreversibleQ
```

### Codon models

```@docs; canonical=false
CodonPartition
```

We offer constructors for the following Q matrix parameterizations:
- `MolecularEvolution.MG94_F3x4` - for example usage, see [`Example 3: FUBAR`](@ref)
- `MolecularEvolution.MG94_F61`
- `MolecularEvolution.HB98_F61`
- `MolecularEvolution.HB98_AAfit`
Use help mode, `?`, in the REPL to find out more about the above.

### Miscellaneous models
```@docs; canonical=false
InterpolatedDiscreteModel
PiQ
```

## Continuous models

The partition of interest is:
```@docs; canonical=false
GaussianPartition
```

And then there are two `BranchModel`s which are compatible with the above partition, namely:
```@docs; canonical=false
BrownianMotion
ZeroDriftOrnsteinUhlenbeck
```

## Compound models

### Branch-wise mixture
```@docs; canonical=false
BWMModel
```

### CAT
```@docs; canonical=false
CATModel
CATPartition
```

### Covarion
```@docs; canonical=false
CovarionModel
CovarionPartition
```

### Site-wise mixture
See [Example 2: GTR+Gamma](@ref) for example usage.
```@docs; canonical=false
SWMModel
SWMPartition
```

## Lazy models

### LazyPartition

```@docs; canonical=false
LazyPartition
```

#### Examples

##### Example 1: Initializing for an upward pass
Now, we show how to wrap the `CodonPartition`s from [Example 3: FUBAR](@ref) with `LazyPartition`:

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

##### Example 2: Initializing for a downward pass
Now, we show how to wrap the `GaussianPartition`s from [Quick example: Likelihood calculations under phylogenetic Brownian motion:](@ref) with `LazyPartition`:

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
!!! note
    Now, we provided a direction for `lazyprep!`. The direction is an instance of `LazyDown`, which was initialized with the `isleafnode` function. The function `isleafnode` dictates if a node saves its sampled observation after a down pass. If you use `direction=LazyDown()`, every node saves its observation.

#### Surrounding LazyPartition
```@docs; canonical=false
lazyprep!
LazyUp
LazyDown
```