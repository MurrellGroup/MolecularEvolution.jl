#This  contains everything that refers to FelNode or to any of the abstract models or partitions.

mutable struct FelNode <: AbstractTreeNode
    parent::Union{FelNode,Nothing}
    children::Array{FelNode,1}
    branchlength::Float64
    name::AbstractString
    nodeindex::Int
    seqindex::Int
    state_path::Vector{StatePath}
    branch_params::Vector{Float64} #Parameters that apply to each branch of the tree - consider nuking this.
    node_data::Dict #For storing anything you like
    message::Vector{Partition} #Vector of Partitions, where each Partition is P(obs|state,model)
    parent_message::Vector{Partition} #This the "downward" message, P(obs,state|model), but at the top of the branch, before it has been propped over the branch.
    child_messages::Vector{Vector{Partition}} #One "up" message" for each child, after prop. Cached to make felsenstein_down! avoid re-computing
    FelNode(branchlength::Float64, name::AbstractString) =
        new(nothing, FelNode[], branchlength, name)
    FelNode() = FelNode(0.0, "")
end

broadcastable(x::FelNode) = Ref(x) #???

function Base.show(io::IO, z::FelNode)
    println(io, "FelNode")
    #println(io, "Type: ", typeof(z.branchlength))
    println(io, "nodeindex: $(z.nodeindex)\nRoot: $(isroot(z))\nLeaf: $(isleafnode(z))")
    println(io, "Defined fields:")
    for f in fieldnames(typeof(z))
        println(io, isdefined(z, f), "\t ", f)
    end
end

function print_traversal(node::FelNode)
    print(node.nodeindex, " ", node.branchlength, " ", node.name)
    if isdefined(node, :branch_params)
        print(" ", node.branch_params)
    end
    if !isleafnode(node)
        println(" ", [n.nodeindex for n in node.children])
        for nod in node.children
            print_traversal(nod)
        end
    else
        println()
    end
end

function set_node_indices!(root::FelNode; starting_index = 1)
    count = 1
    for node in getnodelist(root)
        node.nodeindex = count
        count += 1
    end
end

#Does this need a parametric type to prevent type mismatch?
"""
    internal_message_init!(tree::FelNode, empty_message::Vector{<:Partition})

    Initializes the message template for each node in the tree, allocating space for each partition.
"""
function internal_message_init!(tree::FelNode, empty_message::Vector{<:Partition})
    for node in getnodelist(tree)
        if !isleafnode(node)
            node.child_messages = [deepcopy(empty_message) for i in node.children]
        end
        node.message = deepcopy(empty_message)
        node.parent_message = deepcopy(empty_message)
    end
end

"""
    internal_message_init!(tree::FelNode, partition::Partition)

    Initializes the message template for each node in the tree, as an array of the partition.
"""
function MolecularEvolution.internal_message_init!(tree::FelNode, part_template::Partition)
    internal_message_init!(tree, [part_template])
end

function random_leaf_init!(tree::FelNode, empty_message::Vector{<:Partition})
    for node in getleaflist(tree)
        node.message = deepcopy(empty_message)
        for part in node.message
            for site = 1:part.sites
                part.state[rand(1:part.states), site] = 1.0
            end
        end
    end
end

#Note the different name here...
function mixed_type_equilibrium_message(
    model_vec::Vector{<:BranchModel},
    message_template::Vector{<:Partition},
)
    out_mess = deepcopy(message_template)
    for part = 1:length(message_template)
        out_mess[part] = eq_freq_from_template(model_vec[part], message_template[part])
    end
    return out_mess
end
