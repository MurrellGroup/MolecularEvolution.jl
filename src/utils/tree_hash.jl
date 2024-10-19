#Note to future self/others. The algorthim implemented in this file compares two trees, and outputs a list of nodes
#in the query tree that are not present in the ref tree (in the sense that no node on the ref tree induces the same
#split upon the leaves). If you sort the hash tuples inside get_node_hashes, then the comparison is invariant to 
#re-rooting. If you work through this algorithm and decide "no, it isn't actually invariant to re-rooting - for that 
#to work you'd need to use three hashes per node, one for each branch exiting the node" then you would be very sensible, 
#but also incorrect.

tuple_sort(t::Tuple{UInt64,UInt64}) = ifelse(t[1] < t[2], t, (t[2], t[1]))

function node_hash_split(node, hash_container, node_container, name2hash; push_leaves = false)
    if isleafnode(node)
        if push_leaves
            push!(hash_container, name2hash[node.name])
            push!(node_container, node)
        end
        return name2hash[node.name]
    else
        child_hashes = [
            node_hash_split(nc, hash_container, node_container, name2hash, push_leaves = push_leaves) for
            nc in node.children
        ]
        #merge child hashes with xor as we go up the tree
        first_hash = child_hashes[1]
        for ch in child_hashes[2:end]
            first_hash = xor(first_hash, ch)
        end
        push!(hash_container, first_hash)
        push!(node_container, node)
        return first_hash
    end
end

function get_node_hashes(newt; push_leaves = false)
    leafnames = [n.name for n in getleaflist(newt)]
    leafhashes = hash.(leafnames)
    all_names_hash = xor(leafhashes...)
    name2hash = Dict(zip(leafnames, leafhashes))
    hash_container = UInt64[]
    node_container = FelNode[]
    #This puts things in the containers
    node_hash_split(newt, hash_container, node_container, name2hash, push_leaves = push_leaves)
    #This makes a hash that matches everything except the given node
    other_hash = xor.(hash_container, all_names_hash)
    #Sort these, to make comparisons order invariant, which makes the comparison rooting invariant
    #Consider making this sort an option, and then we can have a rooted comparison and an unrooted one
    sorted_hash_pairs = tuple_sort.(collect(zip(hash_container, other_hash)))
    return sorted_hash_pairs, node_container
end

export tree_diff
#returns nodes in the query that don't have matching splits in the reference
function tree_diff(query, reference)
    newt_hc, newt_nc = get_node_hashes(query)
    n_hc, n_nc = get_node_hashes(reference)
    hashset = Set(n_hc)
    changed_nodes = newt_nc[[!(n in hashset) for n in newt_hc]]
    return changed_nodes
end

export tree_match_pairs
function tree_match_pairs(query, reference; push_leaves = false)
    newt_hc, newt_nc = get_node_hashes(query, push_leaves = push_leaves)
    n_hc, n_nc = get_node_hashes(reference, push_leaves = push_leaves)
    newt_hash2node = Dict(zip(newt_hc, newt_nc))
    n_hash2node = Dict(zip(n_hc, n_nc))
    return map(h -> (newt_hash2node[h], n_hash2node[h]), filter(x -> haskey(n_hash2node, x), newt_hc))
end
#=
function tree_comp(query, reference, condition)
    newt_hc, newt_nc = get_node_hashes(query)
    n_hc, n_nc = get_node_hashes(reference)
    changed_nodes = filter(condition(n_hc), newt_hc)
    return changed_nodes
end

export tree_diff
tree_diff(query, reference) = tree_comp(query, reference, !in ∘ Set)

export tree_match
function tree_match(query, reference)
    filter(newt_hc, n_hc) = map(in(Set(n_hc)), newt_hc)
    return tree_comp(query, reference, filter)
end

tree_match(query, reference) = tree_comp(query, reference, in ∘ Set)
=#