#=
#Move this
export scrape_figtree_colored_nexus
"""
    scrape_figtree_colored_nexus(fname; custom_labels = String[])

Takes a nexus file from FigTree, where branches have been colored. Replaces all color tags with group tags that can be used in the models. Can add custom labels too. Should consider an entire custom dictionary as well in future.
"""
function scrape_figtree_colored_nexus(fname; custom_labels = String[])
    start_substr = "[&R] ";
    lines = readlines(fname);
    treeline = lines[findfirst([occursin(start_substr,l) for l in lines])];
    st = findfirst(start_substr, treeline)[end]
    treestr = treeline[st+1:end]
    R = r"\[\&\!color\=\#\w*\]"
    color_tags = union([string(m.match) for m in eachmatch(R, treestr)])
    d = Dict{String,String}()

    if length(color_tags) == 0
        @warn "No color tags detected."
    else
        if length(custom_labels) == 0
            for i in 1:length(color_tags)
                d[color_tags[i]] = "{G$(i)}"
            end
        else
            for i in 1:length(color_tags)
                d[color_tags[i]] = custom_labels[i]
            end
        end
    end
    for k in keys(d)
        treestr = replace(treestr,k=>d[k])
    end
    return treestr, d
end
=#

function populate_message!(message::Vector{<:Partition}, data)
    if length(message) == 1
        obs2partition!(message[1], data)
    elseif length(message) == length(data) #Might break.
        for i = 1:length(message)
            obs2partition!(message[i], data[i])
        end
    else
        @error "Can't figure out how the data layout corresponds to the message layout. Do it yourself."
    end
end

#This means we need "identity!" to be defined for a partition to use it in the populate_tree! helper functions, if data is missing.
function uninformative_message!(message::Vector{<:Partition})
    identity!.(message)
end

#Sets up a tree, including root and lead messages.
function populate_tree!(
    tree::FelNode,
    starting_message::Vector{<:Partition},
    names,
    data;
    init_all_messages = true,
    tolerate_missing = 1, #0 = error if missing; 1 = warn and set to missing data; 2 = set to missing data
)
    if init_all_messages
        internal_message_init!(tree, starting_message)
    else
        tree.parent_message = deepcopy(starting_message)
    end
    name_dic = Dict(zip(names, 1:length(names)))
    for n in getleaflist(tree)
        if haskey(name_dic, n.name)
            populate_message!(n.message, data[name_dic[n.name]])
        else
            warn_str = n.name * " on tree but not found in names."
            if tolerate_missing == 0
                @error warn_str
            end
            if tolerate_missing == 1
                @warn warn_str
            end
            uninformative_message!(n.message)
        end
    end
end

function populate_tree!(
    tree::FelNode,
    starting_partition::Partition,
    names,
    data;
    init_all_messages = true,
    tolerate_missing = 1,
)
    populate_tree!(
        tree,
        [starting_partition],
        names,
        data,
        init_all_messages = init_all_messages,
        tolerate_missing = tolerate_missing,
    )
end

export expected_subs_per_site
"""
    expected_subs_per_site(Q,mu)

Takes a rate matrix Q and an equilibrium frequency vector, and calculates the expected number of substitutions per site.
"""
function expected_subs_per_site(Q, mu)
    t = 0.0
    mu_norm = mu ./ sum(mu)
    for i = 1:length(mu)
        t += -Q[i, i] * mu_norm[i]
    end
    return t
end

#Consider relocating this.
export tree_polish!
"""
tree_polish!(newt, models; tol = 10^-4, verbose = 1, topology = true)

Takes a tree and a model function, and optimizes branch lengths and, optionally, topology. Returns final LL. Set `verbose=0` to suppress output.
Note: This is not intended for an exhaustive tree search (which requires different heuristics), but rather to polish a tree that is already relatively close to the optimum.
"""
function tree_polish!(newt, models; tol = 10^-4, verbose = 1, topology = true)
    LL = log_likelihood!(newt, models)
    verbose > 0 && println("LL: ", LL)
    d = Inf
    for i = 1:50
        if topology
            felsenstein_down!(newt, models)
            nni_optim!(newt, models)
            log_likelihood!(newt, models)
        end
        felsenstein_down!(newt, models)
        branchlength_optim!(newt, models)
        newLL = log_likelihood!(newt, models)
        verbose > 0 && println("LL: ", newLL)
        d = newLL - LL
        if d < tol
            break
        end
        LL = newLL
    end
    return LL
end


#Replace normalize() in current codebase with this
function sum2one(vec)
    return vec ./ sum(vec)
end

function sum2one!(vec)
    vec .= vec ./ sum(vec)
end

"""
    linear_scale(val,in_min,in_max,out_min,out_max)

Linearly maps val which lives in [in_min,in_max] to a value in [out_min,out_max]
"""
function linear_scale(val, in_min, in_max, out_min, out_max)
    (val - in_min) * ((out_max - out_min) / (in_max - in_min)) + out_min
end

#Draw a categorical value, returned as one-of-k vector.
#We Should likely switch to a 1-hot type, that eg. knows to just index when multiplied onto a matrix.
function one_hot_sample(vec::Vector{<:Real})
    ind = sample(1:length(vec), Weights(vec))
    sampled = zeros(eltype(vec), size(vec))
    sampled[ind] = 1.0
    return sampled
end

function scaled_prob_domain(vec)
    return exp.(vec .- maximum(vec))
end

#The deprecated scale!
function scale_cols_by_vec!(mat, vec)
    rmul!(mat, Diagonal(vec))
end

function name2seq(seqnames::Vector{String}, seqs::Vector{String})
    name2seq_dict = Dict{String,String}()
    for i = 1:length(seqnames)
        name2seq_dict[seqnames[i]] = seqs[i]
    end
    return name2seq_dict
end

function is_one_hot(vec::Vector{T}) where {T<:Real}
    return sum(vec) == maximum(vec) == one(T)
end
function is_one_hot(arr::Array{T,2}) where {T<:Real}
    return all(is_one_hot.(collect.(eachcol(arr))))
end

function binarize!(tree)
    for x in getnodelist(tree)
        while length(x.children) > 2
            new_node = typeof(tree)()
            new_node.children = x.children[1:2]
            x.children = vcat(new_node, x.children[3:end])
            for c in new_node.children
                c.parent = new_node
            end
            new_node.parent = x
        end
    end
end

export read_newick_tree
"""
read_newick_tree(treefile)

Reads in a tree from a file, of type FelNode
"""
function read_newick_tree(
    treefile::String;
    binarize = true,
    ladderize = true,
    strip_single_quotes = true,
)
    treestring = read(treefile, String)
    tree = gettreefromnewick(treestring, FelNode)
    if binarize
        binarize!(tree)
    end
    if ladderize
        ladderize!(tree)
    end
    if strip_single_quotes
        for n in getnodelist(tree)
            n.name = replace(n.name, "'" => "")
        end
    end
    return tree
end
