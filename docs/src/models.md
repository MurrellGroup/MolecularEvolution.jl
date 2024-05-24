# Models

Coming soon.

## Discrete state models

### Codon models

## Continuous models

## Compound models

## Lazy models

### LazyPartition

With this data structure, you can wrap a partition of choice. 
The idea is that in some message passing algorithms, there is only a wave of partitions which need to actualize. 
For instance, a wave following a root-leaf path, or a depth-first traversal.
In which case, we can be more economical with our memory consumption.
With a worst case memory complexity of O(log(n)), where n is the number of nodes, functionality is provided for:
- `log_likelihood!`
- `felsenstein!`
- `sample_down`
(Further functionality is under development)

#### Examples

##### Example 1: Initializing for an upward pass
Now, we show how to wrap the `GappyAminoAcidPartition`s from [Example 1](examples.md#Example-1:-Amino-acid-ancestral-reconstruction-and-visualization) with `LazyPartition`:
You simply swap out these two lines
```julia
initial_partition = GappyAminoAcidPartition(AA_freqs,length(seqs[1]))
populate_tree!(tree,initial_partition,seqnames,seqs)
```

With this
```julia
eq_partition = GappyAminoAcidPartition(AA_freqs,length(seqs[1]))
initial_partition = LazyPartition{GappyAminoAcidPartition}(nothing)
populate_tree!(tree,initial_partition,seqnames,seqs)
lazyprep!(tree, [eq_partition])
```

By this slight modification, we go from initializing and using 212 partitions to 7 during the subsequent `log_likelihood!` calls. There is no significant decrease in performance recorded from this switch.

##### Example 2: Initializing for a downward pass
Now, we show how to wrap the `GaussianPartition`s from [Quick example](../../README.md#quick-example-likelihood-calculations-under-phylogenetic-brownian-motion) with `LazyPartition`:
You simply swap out this line of code
```julia
internal_message_init!(tree, GaussianPartition())
```

With this (technically we only add 1 LOC)
```julia
initial_partition = GaussianPartition()
lazy_initial_partition = LazyPartition{GaussianPartition}(nothing)
internal_message_init!(tree, lazy_initial_partition)
lazyprep!(tree, initial_partition, direction=LazyDown(isleafnode))
```

Note that we now provided a direction for `lazyprep!`. The direction is an instance of `LazyDown`, which was initialized with the `isleafnode` function. The function `isleafnode` dictates if a node saves its sampled observation after a down pass. If you use `direction=LazyDown()`, every node saves its observation.