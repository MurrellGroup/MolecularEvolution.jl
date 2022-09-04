#Note to future self/others. The algorthim implemented in this file compares two trees, and outputs a list of nodes
#in the query tree that are not present in the ref tree (in the sense that no node on the ref tree induces the same
#split upon the leaves). If you sort the hash tuples inside get_node_hashes, then the comparison is invariant to 
#re-rooting. If you work through this algorithm and decide "no, it isn't actually invariant to re-rooting - for that 
#to work you'd need to use three hashes per node, one for each branch exiting the node" then you would be very sensible, 
#but also incorrect.

tuple_sort(t::Tuple{UInt64, UInt64}) = ifelse(t[1] < t[2], t, (t[2],t[1]))

function node_hash_split(node, hash_container,node_container, name2hash)
    if isleafnode(node)
        #consider if we want the children pairs too, and push in here too
        return name2hash[node.name]
    else
    child_hashes = [node_hash_split(nc,hash_container,node_container, name2hash) for nc in node.children]
    #merge child hashes with xor as we go up the tree
    first_hash = child_hashes[1]
    for ch in child_hashes[2:end]
        first_hash = xor(first_hash,ch)
    end
    push!(hash_container,first_hash)
    push!(node_container,node)
    return first_hash
    end
end

function get_node_hashes(newt)
    leafnames = [n.name for n in getleaflist(newt)];
    leafhashes = hash.(leafnames)
    all_names_hash = xor(leafhashes...);
    name2hash = Dict(zip(leafnames,leafhashes));
    hash_container = UInt64[]
    node_container = FelNode[]
    #This puts things in the containers
    node_hash_split(newt, hash_container, node_container, name2hash)
    #This makes a hash that matches everything except the given node
    other_hash = xor.(hash_container,all_names_hash)
    #Sort these, to make comparisons order invariant, which makes the comparison rooting invariant
    #Consider making this sort an option, and then we can have a rooted comparison and an unrooted one
    sorted_hash_pairs = tuple_sort.(collect(zip(hash_container,other_hash)))
    return sorted_hash_pairs, node_container
end

export tree_diff
#returns nodes in the query that don't have matching splits in the reference
function tree_diff(query, reference)
    newt_hc,newt_nc = get_node_hashes(query);
    n_hc,n_nc = get_node_hashes(reference)
    hashset = Set(n_hc);
    changed_nodes = newt_nc[[!(n in hashset) for n in newt_hc]];
    return changed_nodes
end
