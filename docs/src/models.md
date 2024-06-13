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