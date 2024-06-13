# Models

Coming soon.

## Discrete state models

### Codon models

## Continuous models

## Compound models

## Lazy models

### LazyPartition

```@docs; canonical=false
LazyPartition
```

#### Examples

##### Example 1: Initializing for an upward pass
Now, we show how to wrap the `GappyAminoAcidPartition`s from [Example 1: Amino acid ancestral reconstruction and visualization](@ref) with `LazyPartition`:

You simply swap out these two lines
```julia
initial_partition = GappyAminoAcidPartition(AA_freqs,length(seqs[1]))
populate_tree!(tree,initial_partition,seqnames,seqs)
```

With this
```julia
eq_partition = GappyAminoAcidPartition(AA_freqs,length(seqs[1]))
initial_partition = LazyPartition{GappyAminoAcidPartition}()
populate_tree!(tree,initial_partition,seqnames,seqs)
lazyprep!(tree, [eq_partition])
```

By this slight modification, we go from initializing and using 212 partitions to 7 during the subsequent `log_likelihood!` calls. There is no significant decrease in performance recorded from this switch.

##### Example 2: Initializing for a downward pass
Now, we show how to wrap the `GaussianPartition`s from [Quick example: Likelihood calculations under phylogenetic Brownian motion:](@ref) with `LazyPartition`:

You simply swap out this line of code
```julia
internal_message_init!(tree, GaussianPartition())
```

With this (technically we only add 1 LOC)
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