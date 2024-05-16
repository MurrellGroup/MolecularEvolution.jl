# Models

Coming soon.

## Discrete state models

### Codon models

## Continuous models

## Compound models

## Lazy models

### LazyPartition

With this data structure, you can wrap a partition of choice. With a worst case memory complexity of O(log(n)), where n is the number of nodes, functionality is provided for:
- `log_likelihood!`
- `felsenstein!`
(Further functionality is under development)

#### Example
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