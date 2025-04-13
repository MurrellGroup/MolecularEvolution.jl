export HIPSTR
"""
    HIPSTR(trees::Vector{FelNode}; set_branchlengths = true)

Construct a Highest Independent Posterior Subtree Reconstruction (HIPSTR) tree
from a collection of trees.

Returns a single FelNode representing the HIPSTR consensus tree.

If `set_branchlengths = true`, the branch length of a node in the HIPSTR tree will be set to the mean branch length of all nodes from the input trees that have the same clade. (By the same clade, we mean that the set of leaves below the node is the same.) Otherwise, the root branch length is 0.0 and the rest 1.0.

Source: https://www.biorxiv.org/content/10.1101/2024.12.08.627395v1.full.pdf
"""
function HIPSTR(trees::Vector{FelNode}; set_branchlengths = true)
    @info "Starting HIPSTR construction from $(length(trees)) trees"
    
    # Step 1: Collect all clades, their frequencies, and child pairs
    clades_stats = Dict{Tuple{UInt64, UInt64}, CladeStats}()
    leaf_names = Set{String}()
    
    # First identify all leaf nodes across all trees
    for tree in trees
        for leaf in getleaflist(tree)
            push!(leaf_names, leaf.name)
        end
    end
    
    # Create a mapping from leaf names to indices for consistent hashing
    leaf_dict = Dict(name => i for (i, name) in enumerate(sort(collect(leaf_names))))
    
    # Process each tree to collect clade information
    for tree in trees
        collect_clades!(tree, clades_stats, leaf_dict)
    end
    
    # Scale clade frequencies to get credibilities
    n_trees = length(trees)
    for (_, stats) in clades_stats
        stats.frequency /= n_trees
    end
    
    # Step 2: Compute the root clade hash (all tips)
    all_tips = BitSet(1:length(leaf_dict))
    root_hash = hash_clade(all_tips)
    @info "Root clade hash: $root_hash"
    
    # Step 3: Build the credibility cache through post-order traversal
    cred_cache = Dict{Tuple{UInt64, UInt64}, Tuple{Float64, Tuple{UInt64, UInt64}, Tuple{UInt64, UInt64}}}()
    compute_credibility = function(clade_hash)
        # Return from cache if available
        haskey(cred_cache, clade_hash) && return cred_cache[clade_hash][1]
        
        # Base case: single tip or clade not found
        if !haskey(clades_stats, clade_hash) || isempty(clades_stats[clade_hash].child_pairs)
            cred_cache[clade_hash] = (clades_stats[clade_hash].frequency, (0, 0), (0, 0))
            return clades_stats[clade_hash].frequency
        end
        
        # Find the best child pair
        best_cred = 0.0
        best_left = (0, 0)
        best_right = (0, 0)
        
        for (left_hash, right_hash) in clades_stats[clade_hash].child_pairs
            left_cred = compute_credibility(left_hash)
            right_cred = compute_credibility(right_hash)
            
            # Product of the credibilities
            pair_cred = left_cred * right_cred * clades_stats[clade_hash].frequency
            
            if pair_cred > best_cred
                best_cred = pair_cred
                best_left = left_hash
                best_right = right_hash
            end
        end
        
        # Cache and return
        cred_cache[clade_hash] = (best_cred, best_left, best_right)
        return best_cred
    end
    
    # Compute credibility for the root clade
    compute_credibility(root_hash)
    
    # Step 4: Construct the HIPSTR tree through another traversal
    reverse_leaf_dict = Dict(i => name for (name, i) in leaf_dict)
    
    # Function to build the tree recursively
    function build_tree(clade_hash)
        _, left_hash, right_hash = cred_cache[clade_hash]
        
        # Handle leaf case
        if left_hash == (0, 0) && right_hash == (0, 0)
            # Determine which tip this is
            for (index, name) in reverse_leaf_dict
                tip_hash = hash_clade(BitSet([index]))
                if tip_hash == clade_hash
                    node = FelNode(1.0, name)
                    node.seqindex = index
                    return node
                end
            end
            error("Failed to find leaf corresponding to hash $clade_hash")
        end
        
        # Internal node
        node = FelNode(1.0, "")
        
        # Add children
        left_child = build_tree(left_hash)
        right_child = build_tree(right_hash)
        
        # Default branch lengths to 1.0 if we don't have better information
        left_child.branchlength = 1.0
        right_child.branchlength = 1.0
        
        left_child.parent = node
        right_child.parent = node
        push!(node.children, left_child)
        push!(node.children, right_child)
        
        return node
    end
    
    # Build the final tree
    hipstr_tree = build_tree(root_hash)
    # Set the root branchlength to 0.0
    hipstr_tree.branchlength = 0.0
    
    # Set node indices
    set_node_indices!(hipstr_tree)

    # Set branch lengths to mean
    set_branchlengths && set_mean_branchlengths!(hipstr_tree, trees)
    
    return hipstr_tree
end

function set_mean_branchlengths!(tree::FelNode, trees::Vector{FelNode})
    #Set branchlengths to 0.0, we'll use these as accumulators and then eventually normalize
    tree_nodes = nodes(tree)
    for node in tree_nodes
        node.branchlength = 0.0
    end
    #Initialize container for counting the amount of terms in the accumulators
    branch_length_counts = Dict(zip(tree_nodes, zeros(Int64, length(tree_nodes))))
    #Go through all the trees, t...
    for t in trees
        #and for each matching clade, (tree_node, t_node), between tree and t...
        matching_pairs = tree_match_pairs(tree, t, push_leaves = true)
        for (tree_node, t_node) in matching_pairs
            #add the branch length of t_node in t to the branch length of tree_node in tree...
            tree_node.branchlength += t_node.branchlength
            #and increment the counter for tree_node.
            branch_length_counts[tree_node] += 1
        end
    end
    #Then normalize the branch length of all nodes in tree by the number of matching clades.
    for node in tree_nodes
        branch_length_counts[node] == 0 && @warn "Branch length counts for node $(node.name) with nodeindex $(node.nodeindex) is 0. Coming from a `HIPSTR` call, this number should be strictly positive."
        node.branchlength /= max(1, branch_length_counts[node])
    end
    #^^ This deviates somewhat from the HIPSTR paper. They find distributions over node ages (with respect to a molecular clock), which differs a bit from the notion of branch length.
end
"""
Store statistics about a clade: its frequency and observed child pairs.
"""
mutable struct CladeStats
    frequency::Float64
    child_pairs::Set{Tuple{Tuple{UInt64, UInt64}, Tuple{UInt64, UInt64}}}
    
    CladeStats() = new(0.0, Set{Tuple{Tuple{UInt64, UInt64}, Tuple{UInt64, UInt64}}}())
end

"""
Compute a hash for a clade based on its tips.
"""
function hash_clade(tips::BitSet)
    h1 = hash(tips)
    h2 = hash(reverse(collect(tips)))
    return (h1, h2)
end

"""
Recursively collect clades from a tree.
"""
function collect_clades!(node::FelNode, clades_stats::Dict{Tuple{UInt64, UInt64}, CladeStats}, leaf_dict::Dict{String, Int})
    # Get tips under this node
    tips = BitSet()
    
    if isleafnode(node)
        # For a leaf, the tips are just this node
        if haskey(leaf_dict, node.name)
            push!(tips, leaf_dict[node.name])
        else
            # Skip if the leaf name is not recognized
            @warn "Skipping unrecognized leaf name: $(node.name)"
            return tips
        end
    else
        # For internal nodes, combine tips from children
        for child in node.children
            union!(tips, collect_clades!(child, clades_stats, leaf_dict))
        end
    end
    
    # Compute the clade hash
    clade_hash = hash_clade(tips)
    
    # Update clade stats
    if !haskey(clades_stats, clade_hash)
        clades_stats[clade_hash] = CladeStats()
    end
    clades_stats[clade_hash].frequency += 1
    
    # For internal nodes, record child pairs
    if !isleafnode(node) && length(node.children) == 2
        left_tips = BitSet()
        for leaf in getleaflist(node.children[1])
            if haskey(leaf_dict, leaf.name)
                push!(left_tips, leaf_dict[leaf.name])
            end
        end
        
        right_tips = BitSet()
        for leaf in getleaflist(node.children[2])
            if haskey(leaf_dict, leaf.name)
                push!(right_tips, leaf_dict[leaf.name])
            end
        end
        
        left_hash = hash_clade(left_tips)
        right_hash = hash_clade(right_tips)
        
        # Add the child pair to the set for this clade
        push!(clades_stats[clade_hash].child_pairs, (left_hash, right_hash))
    end
    
    return tips
end
