"Internal function. Helper for bfs_mapreduce and dfs_mapreduce"
function _mapreduce(
    start_node::AbstractTreeNode,
    map_reduce::T,
    aggregator,
    queue_ordering_func!,
) where {T<:Function}
    queue = [(start_node, start_node)]

    while length(queue) > 0
        curr_node, prev_node = queue_ordering_func!(queue)
        parent_and_children =
            filter(x -> !isnothing(x), vcat(curr_node.parent, curr_node.children))
        new_queue_elems = [
            (next_node, curr_node) for
            next_node in filter(x -> x != prev_node, parent_and_children)
        ]
        #println(queue)
        #println(new_queue_elems)
        append!(queue, new_queue_elems)
        map_reduce(curr_node, prev_node, aggregator)
    end
end


export bfs_mapreduce
"""Performs a BFS map-reduce over the tree, starting at a given node
For each node, map_reduce is called as: 
  map_reduce(curr_node::FelNode, prev_node::FelNode, aggregator)
where prev_node is the previous node visited on the path from the start node to the current node
It is expected to update the aggregator, and not return anything.

Not exactly conventional map-reduce, as map-reduce calls may rely on state in the aggregator added
by map-reduce calls on other nodes visited earlier.
"""
function bfs_mapreduce(
    start_node::AbstractTreeNode,
    map_reduce::T,
    aggregator,
) where {T<:Function}
    _mapreduce(start_node, map_reduce, aggregator, Base.popfirst!)
end


export dfs_mapreduce
"""Performs a DFS map-reduce over the tree, starting at a given node
See bfs_mapreduce for more details.
"""
function dfs_mapreduce(
    start_node::AbstractTreeNode,
    map_reduce::T,
    aggregator,
) where {T<:Function}
    _mapreduce(start_node, map_reduce, aggregator, Base.pop!)
end

export node_distances
"Compute the distance to all other nodes from a given node"
function node_distances(from_node::AbstractTreeNode)::Dict{<:AbstractTreeNode,Float64}
    function map_reduce(
        curr_node::AbstractTreeNode,
        prev_node::AbstractTreeNode,
        aggregator::Dict{<:AbstractTreeNode,Float64},
    )
        bl = curr_node.parent == prev_node ? curr_node.branchlength : prev_node.branchlength
        aggregator[curr_node] = aggregator[prev_node] + bl
    end
    aggregator = Dict{FelNode,Float64}(from_node => 0)
    bfs_mapreduce(from_node, map_reduce, aggregator)
    return aggregator
end

"Provides a list of parent nodes nodes from this node up to the root node"
function parent_list(node::FelNode)::Array{<:FelNode,1}
    list = [node]
    while !isnothing(node.parent)
        push!(list, node.parent)
        node = node.parent
    end
    return list
end

export shortest_path_between_nodes
"""
Shortest path between nodes, returned as two lists, each starting with one of the two nodes, 
and ending with the common ancestor
"""
function shortest_path_between_nodes(
    node1::FelNode,
    node2::FelNode,
)::Tuple{Array{<:FelNode,1},Array{<:FelNode,1}}
    parent_list1 = parent_list(node1)
    parent_list2 = parent_list(node2)

    i = 0
    while i < length(parent_list1) &&
              i < length(parent_list2) &&
              parent_list1[end-i] == parent_list2[end-i]
        i += 1
    end

    return parent_list1[1:end-(i-1)], parent_list2[1:end-(i-1)]
end

export longest_path
"""
Returns the longest path in a tree
For convenience, this is returned as two lists of form:
    [leaf_node, parent_node, .... root]
Where the leaf_node nodes are selected to be the furthest away
"""
function longest_path(root::FelNode)
    dists = node_distances(root)
    #Pick the maximum distance, breaking ties using rand()
    #This is needed to prevent a comparison between nodes which is undefined
    _, _, furthest_node = maximum(x -> (x[2], rand(), x[1]), dists)

    dists_from_furthest_node = node_distances(furthest_node)

    _, _, other_furthest_node = maximum(x -> (x[2], rand(), x[1]), dists_from_furthest_node)

    return shortest_path_between_nodes(furthest_node, other_furthest_node)
end

export midpoint
"Returns a midpoint as a node and a distance above it where the midpoint is"
function midpoint(root::FelNode)::Tuple{<:FelNode,Float64}
    path1, path2 = longest_path(root)

    #compute the total branch lengths down from the root to either farthest node
    side_length1 = sum(x -> x.branchlength, path1[1:end-1])
    side_length2 = sum(x -> x.branchlength, path2[1:end-1])

    total_length = side_length1 + side_length2

    farther_side = side_length1 > side_length2 ? path1 : path2
    remaining_length =
        (total_length / 2) - (side_length1 > side_length2 ? side_length2 : side_length1)

    for node in Iterators.reverse(farther_side[1:end-1])
        remaining_length -= node.branchlength
        if remaining_length < 0
            return (node, -remaining_length)
        end
    end
end
