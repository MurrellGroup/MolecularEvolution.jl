
# Input/Output {#Input/Output}
<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.write_nexus' href='#MolecularEvolution.write_nexus'><span class="jlbinding">MolecularEvolution.write_nexus</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
write_nexus(fname::String,tree::FelNode)
```


Writes the tree as a nexus file, suitable for opening in eg. FigTree. Data in the `node_data` dictionary will be converted into annotations. Only tested for simple `node_data` formats and types.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/utils/misc.jl#L281-L287" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.newick' href='#MolecularEvolution.newick'><span class="jlbinding">MolecularEvolution.newick</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
newick(root)
```


Returns a newick string representation of the tree.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/core/nodes/AbstractTreeNode.jl#L626-L630" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.read_newick_tree' href='#MolecularEvolution.read_newick_tree'><span class="jlbinding">MolecularEvolution.read_newick_tree</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



read_newick_tree(treefile)

Reads in a tree from a file, of type FelNode


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/utils/misc.jl#L256-L260" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.populate_tree!' href='#MolecularEvolution.populate_tree!'><span class="jlbinding">MolecularEvolution.populate_tree!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
populate_tree!(tree::FelNode, starting_message, names, data; init_all_messages = true, tolerate_missing = 1, leaf_name_transform = x -> x)
```


Takes a tree, and a `starting_message` (which will serve as the memory template for populating messages all over the tree). `starting_message` can be a message (ie. a vector of Partitions), but will also work with a single Partition (although the tree) will still be populated with a length-1 vector of Partitions. Further, as long as `obs2partition` is implemented for your Partition type, the leaf nodes will be populated with the data from `data`, matching the names on each leaf. When a leaf on the tree has a name that doesn&#39;t match anything in `names`, then if
- `tolerate_missing = 0`, an error will be thrown
  
- `tolerate_missing = 1`, a warning will be thrown, and the message will be set to the uninformative message (requires identity!(::Partition) to be defined)
  
- `tolerate_missing = 2`, the message will be set to the uninformative message, without warnings (requires identity!(::Partition) to be defined)
  

A renaming function that can eg. strip tags from the tree when matching leaf names with `names` can be passed to `leaf_name_transform`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/utils/misc.jl#L110-L122" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.read_fasta' href='#MolecularEvolution.read_fasta'><span class="jlbinding">MolecularEvolution.read_fasta</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
read_fasta(filepath::String)
```


Reads in a fasta file and returns a tuple of (seqnames, seqs).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/utils/fasta_io.jl#L4-L8" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MolecularEvolution.write_fasta' href='#MolecularEvolution.write_fasta'><span class="jlbinding">MolecularEvolution.write_fasta</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
write_fasta(filepath::String, sequences::Vector{String}; seq_names = nothing)
```


Writes a fasta file from a vector of sequences, with optional seq_names.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MurrellGroup/MolecularEvolution.jl/blob/db088346584c687f47f57a4fc109427cadc63e91/src/utils/fasta_io.jl#L18-L22" target="_blank" rel="noreferrer">source</a></Badge>

</details>

