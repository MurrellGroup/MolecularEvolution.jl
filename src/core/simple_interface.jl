#Making this broadcast if called on a message instead of a partition.
partition2obs(vp::Vector{<:Partition}) = partition2obs.(vp)

"""
    leaf_samples(tree; partition_inds = 1)

Returns the result of `partition2obs` for each leaf of the tree. Can be used eg. after `sample_down!` is called.
If using a eg. codon model, this will extract a string from the CodonPartition on each leaf.
Acts upon the first partition by default, but this can be changed by setting `partition_inds`, which can also be a vector of indices, 
in which case the result will be a vector for each leaf.
"""
leaf_samples(tree::FelNode; partition_inds = 1) = [partition2obs(n.message[partition_inds]) for n in getleaflist(tree)]

"""
    node_samples(tree; partition_inds = 1)

Returns the result of `partition2obs` for each node of the tree (including internal nodes, and the root). Can be used eg. after `sample_down!` is called.
If using a eg. codon model, this will extract a string from the CodonPartition on each node.
Acts upon the first partition by default, but this can be changed by setting `partition_inds`, which can also be a vector of indices, 
in which case the result will be a vector for each node.
"""
node_samples(tree::FelNode; partition_inds = 1) = [partition2obs(n.message[partition_inds]) for n in getnodelist(tree)]

"""
    leaf_names(tree)

Returns the names of the leaves of the tree.
"""
leaf_names(tree::FelNode) = [n.name for n in getleaflist(tree)]

"""
    node_names(tree)

Returns the names of the nodes of the tree.
"""
node_names(tree::FelNode) = [n.name for n in getnodelist(tree)]

"""
    allocate!(tree, partition_or_message)

Allocates initial messages for all nodes in the tree, copying the passed-in message template.
If passed a partition, then this will assume the message template is a vector containing just that partition.
"""    
allocate!(tree, m_or_p) = internal_message_init!(tree, m_or_p)

"""
    leaves(tree)

Returns the leaves of the tree, as a vector.
"""
leaves(tree) = getleaflist(tree)

"""
    nodes(tree)

Returns the nodes of the tree (including internal nodes and the root), as a vector.
"""
nodes(tree) = getnodelist(tree)

"""
    internal_nodes(tree)

Returns the internal nodes of the tree (including the root), as a vector.
"""
internal_nodes(tree) = getnonleaflist(tree)


