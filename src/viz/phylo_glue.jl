#If you have Phylo and Plots installed, you can use this to draw trees with their code.

#Need to add a version that migrates the node_data Dict to the Phylo representation
#data_function must return a key,value pair
function add_node_to_phylo_tree(
    molev_node::FelNode,
    phylo_tree,
    phylo_node,
    node_counter;
    copy_dict = true,
    data_function = (x -> Tuple{String,Float64}[]),
)
    for tup in data_function(molev_node)
        key, value = tup
        phylo_node.data[key] = value
    end

    if isdefined(molev_node, :node_data)
        for k in keys(molev_node.node_data)
            phylo_node.data[k] = molev_node.node_data[k]
        end
    end

    phylo_node_name = Phylo.getnodename(phylo_tree, phylo_node)
    for c in molev_node.children
        if c.name != ""
            new_child = Phylo.createnode!(phylo_tree, c.name)
        else
            new_child = Phylo.createnode!(phylo_tree)
        end
        Phylo.createbranch!(phylo_tree, phylo_node_name, new_child, c.branchlength)
        add_node_to_phylo_tree(
            c,
            phylo_tree,
            new_child,
            node_counter,
            data_function = data_function,
        )
    end
end

export get_phylo_tree

"""
    get_phylo_tree(molev_root::FelNode; data_function = (x -> Tuple{String,Float64}[]))

Converts a FelNode tree to a Phylo tree. The `data_function` should return a list of tuples of the form (key, value) to be added to the Phylo tree `data` Dictionary.
Any key/value pairs on the FelNode `node_data` Dict will also be added to the Phylo tree.
"""
function get_phylo_tree(molev_root::FelNode; data_function = (x -> Tuple{String,Float64}[]))
    phylo_tree = Phylo.RootedTree()
    node_counter = 1
    if molev_root.name != ""
        root_name = molev_root.name
    else
        root_name = "root"
    end
    phylo_root = Phylo.createnode!(phylo_tree, root_name)
    add_node_to_phylo_tree(
        molev_root,
        phylo_tree,
        phylo_root,
        1,
        data_function = data_function,
    )
    return phylo_tree
end

export values_from_phylo_tree

"""
    values_from_phylo_tree(phylo_tree, key)

    Returns a list of values from the given key in the nodes of the phylo_tree, in an order that is somehow compatible with the order the nodes get plotted in.

"""
function values_from_phylo_tree(phylo_tree, key)
    return [Phylo.getnodedata(phylo_tree, n)[key] for n in Phylo.getnodes(phylo_tree)]
end

export savefig_tweakSVG

"""
    savefig_tweakSVG(fname, plot<:Plots.Plot; hack_bounding_box = true, new_viewbox = nothing, linecap_round = true)

Note: Might only work if you're using the GR backend!!
Saves a figure created using the `Phylo` `Plots` recipe, but tweaks the SVG after export.
`new_viewbox` needs to be an array of 4 numbers, typically starting at `[0 0 plot_width*4 plot_height*4]`
but this lets you add shifts, in case the plot is getting cut off.

eg. `savefig_tweakSVG("export.svg",pl, new_viewbox = [-100, -100, 3000, 4500])`
"""
function savefig_tweakSVG(
    fname, plot::Plots.Plot;
    hack_bounding_box = true,
    new_viewbox = nothing,
    linecap_round = true,
)
    if !(fname[end-3:end] == ".svg")
        fname = fname * ".svg"
    end
    Plots.savefig(plot, fname)
    s = read(fname, String)
    if hack_bounding_box
        s = replace(s, r"<clipPath.*\n.*\n.*</clipPath>" => "")
    end
    if !isnothing(new_viewbox)
        vb = new_viewbox
        s = replace(s, r"viewBox=.*\"" => "viewBox=\"$(vb[1]) $(vb[2]) $(vb[3]) $(vb[4])\"")
    end
    if linecap_round
        s = replace(s, "stroke-linecap:butt" => "stroke-linecap:round")
    end
    open(fname, "w") do f
        write(f, s)
    end
end
