struct InternalPlotAttributes
    y_coords::Dict{String, <:Real}
    node_size::Real
    font_size::Real
    margin::Real
end

function getdistfromrootdict(node::FelNode, xdict::Dict{FelNode, Float64} = Dict{FelNode, Float64}())
    x = node.branchlength
    if isroot(node)
        xdict[node] = x
    else
        xdict[node] = x + xdict[node.parent]
    end
    for child in node.children
        getdistfromrootdict(child, xdict)
    end
    return xdict
end

function WLS(tree, inf_tree, weight_fn; opt_scale = true)
    xdict, inf_xdict = getdistfromrootdict(tree), getdistfromrootdict(inf_tree)
    matched_pairs = tree_match_pairs(tree, inf_tree, push_leaves = true)
    #Make sure root is added to the matched_pairs, but not twice...
    all(x -> !isroot(x[1]), matched_pairs) && push!(matched_pairs, (tree, inf_tree))
    n = length(matched_pairs)
    w, x, y = zeros(n), zeros(n), zeros(n)
    for (i, (node, inf_node)) in enumerate(matched_pairs)
        @assert (w[i] = weight_fn(node)) >= 0
        x[i] = xdict[node]
        y[i] = inf_xdict[inf_node]
    end
    if opt_scale
        sqrt_w = sqrt.(w)
        location, scale = (hcat(sqrt_w, sqrt_w .* x)) \ (sqrt_w .* y)
        if scale > 0.0
            return location, scale
        end
        #If the optimal scale is not well-behaved, we fall back to the opt_scale = false problem
    end
    location, scale = sum(w .* (y - x)) / sum(w), 1.0
    return location, scale
end

# Function to recursively plot internal nodes and edges
function plot_internal!(node, x_parent, scale, line_width, line_alpha, flag, attr)
    if isleafnode(node)
        y = attr.y_coords[node.name]
        x = x_parent + node.branchlength * scale
        if flag
            Plots.scatter!([x], [y], markersize=attr.node_size, color=:blue)
            Plots.annotate!(x + attr.margin, y, Plots.text(node.name, attr.font_size, :left))
        end
        return x, y
    else
        node_x = x_parent + node.branchlength * scale
        children = node.children
        child_coords = [plot_internal!(child, node_x, scale, line_width, line_alpha, flag, attr) for child in children]
        node_y = sum(last.(child_coords)) / length(child_coords)

        for (child_x, child_y) in child_coords
            Plots.plot!([node_x, child_x], [node_y, child_y],
            color=Plots.RGBA(0, 0, 0, line_alpha),
            linewidth=line_width)
        end

        return node_x, node_y
    end
end

export plot_multiple_trees
"""
    plot_multiple_trees(trees, inf_tree; <keyword arguments>)

Plots multiple phylogenetic trees against a reference tree, `inf_tree`.
For each **tree** in `trees`, a linear Weighted Least Squares (WLS) problem (parameterized by the `weight_fn` keyword) is solved for the x-positions of the matching nodes between `inf_tree` and **tree**.

# Keyword Arguments
- `node_size=4`: the size of the nodes in the plot.
- `line_width=0.5`: the width of the branches from `trees`.
- `font_size=10`: the font size for the leaf labels.
- `margin=1.5`: the margin between a leaf node and its label.
- `line_alpha=0.05`: the transparency level of the branches from `trees`.
- `weight_fn=n::FelNode -> ifelse(isroot(n), 1.0, 0.0))`: a function that assigns a weight to a node for the WLS problem.
- `opt_scale=true`: whether to include a scaling parameter for the WLS problem.
"""
function plot_multiple_trees(trees, inf_tree;
    node_size=4,
    line_width=0.5,
    font_size=10,
    margin=1.5,
    line_alpha=0.05,
    weight_fn=n::FelNode -> ifelse(isroot(n), 1.0, 0.0),
    opt_scale=true)
    # Assume all trees have the same leaf set
 
    leaves = getleaflist(inf_tree)
    n_leaves = length(leaves)
    target_max_depth = 100.0
 
    # Initialize plot
    p = Plots.plot(size=(800, 600), legend=false, grid=false, ticks=false, border=:none)
 
    # Calculate y-coordinates for leaves (use the same order for all trees)
    y_coords = Dict(leaf.name => i for (i, leaf) in enumerate(reverse(leaves)))
 
 
    # Find the maximum x-coordinate (tree depth) across all trees
    target_scale = target_max_depth / maximum(values(getdistfromrootdict(inf_tree)))

    #Init global internal plot attributes
    attr = InternalPlotAttributes(y_coords, node_size, font_size, margin / target_scale)
 
    # Plot all trees
    for tree in trees
        ladderize!(tree)
        location, scale = WLS(tree, inf_tree, weight_fn; opt_scale = opt_scale)
        plot_internal!(tree, location, scale, line_width, line_alpha, false, attr)
    end
 
    ladderize!(inf_tree)
    plot_internal!(inf_tree, 0.0, 1.0, 2.0, 1.0, true, attr)
 
    Plots.ylims!(0, n_leaves + 1)
 
    return p
end