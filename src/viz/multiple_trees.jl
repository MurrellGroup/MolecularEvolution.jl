totalbl(tree) = sum([n.branchlength for n in getnodelist(tree)])

struct InternalPlotAttributes
    y_coords::Dict{String, <:Real}
    node_size::Real
    font_size::Real
    margin::Real
end

#I assume equal topology
function getb(inf_node, node, x_inf_parent, x_parent)
    x_inf = x_inf_parent + inf_node.branchlength
    x = x_parent + node.branchlength
    return mapreduce(nodes -> getb(nodes..., x_inf, x), +, zip(inf_node.children, node.children); init=0.0) + 2 * (x - x_inf)
end

function least_squares_δ(inf_tree, tree, x_inf, x)
    a = length(getnodelist(inf_tree))
    b = getb(inf_tree, tree, x_inf, x)
    return -b / (2*a)
end

# Function to recursively plot internal nodes and edges
function plot_internal!(node, x_parent, line_width, line_alpha, flag, attr)
    if isleafnode(node)
        y = attr.y_coords[node.name]
        x = x_parent + node.branchlength
        if flag
            Plots.scatter!([x], [y], markersize=attr.node_size, color=:blue)
            Plots.annotate!(x + attr.margin, y, Plots.text(node.name, attr.font_size, :left))
        end
        return x, y
    else
        children = node.children
        child_coords = [plot_internal!(child, x_parent + node.branchlength, line_width, line_alpha, flag, attr) for child in children]
        node_x = x_parent + node.branchlength
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
function plot_multiple_trees(trees, inf_tree;
    node_size=4,
    line_width=0.5,
    font_size=10,
    margin=1.0,
    line_alpha=0.05,
    least_squares=false)
    # Assume all trees have the same leaf set
 
    leaves = getleaflist(inf_tree)
    n_leaves = length(leaves)
    target_total_bl = 100.0
 
    # Initialize plot
    p = Plots.plot(size=(800, 600), legend=false, grid=false, ticks=false, border=:none)
 
    # Calculate y-coordinates for leaves (use the same order for all trees)
    y_coords = Dict(leaf.name => i for (i, leaf) in enumerate(reverse(leaves)))
 
 
    # Find the maximum x-coordinate (tree depth) across all trees
    tbl = totalbl(inf_tree)
    for n in getnodelist(inf_tree)
        n.branchlength *= target_total_bl/tbl
    end
    #max_depth = maximum(values(node_distances(inf_tree)))

    #Init global internal plot attributes
    attr = InternalPlotAttributes(y_coords, node_size, font_size, margin)
 
    # Plot all trees
    for tree in trees
        tbl = totalbl(tree)
        for n in getnodelist(tree)
            n.branchlength *= target_total_bl/tbl
        end
        ladderize!(tree)
        δ = least_squares ? least_squares_δ(inf_tree, tree, 0.0, 0.0) : 0.0
        plot_internal!(tree, δ, line_width, line_alpha, false, attr)
    end
 
    ladderize!(inf_tree)
    plot_internal!(inf_tree, 0, 2.0, 1.0, true, attr)
 
    Plots.ylims!(0, n_leaves + 1)
 
    return p
end