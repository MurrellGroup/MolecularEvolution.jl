#=Required fields in inheritors:
  parent::Nullable{TreeNode}
  children::Array{TreeNode,1}
  branchlength::Float64
  name::AbstractString
  seqindex::Int

  Required constructors in inheritors:
  T()
=#

export AbstractTreeNode
abstract type AbstractTreeNode end

Base.eltype(node::T) where T <: AbstractTreeNode = T

#Currently disagree with these:
#=
Base.length(node::T) where T <: AbstractTreeNode = length(node.children)
Base.iterate(node::T, state=1) where T <: AbstractTreeNode = return state > length(node.children) ? nothing : node.children[state], state+1
=#

function print_traversal(node::AbstractTreeNode)
    print(node.nodeindex," ", node.branchlength, " ", node.name)
    if !isleafnode(node)
        println(" ",[n.nodeindex for n in node.children])
        for nod in node.children
            print_traversal(nod)
        end
    else
        println()
    end
end

export roottree
function roottree(root::T,index::Int=1) where T <: AbstractTreeNode
  #Call constructor of node type
  newroot = T()
  child = splice!(root.children,index)
  child.parent = Union{T,Nothing}()
  addchild(newroot,child)
  br = child.branchlength/2.0
  child.branchlength = br
  root.branchlength = br
  addchild(newroot,root)
  return newroot
end

export isroot
function isroot(node::T) where T <: AbstractTreeNode
  return node.parent === nothing
end

export isleafnode
function isleafnode(node::T) where T <: AbstractTreeNode
  return length(node.children) == 0
end


export istreeconsistent
"""
	istreeconsistent(root)

Checks whether the `:parent` field is set to be consistent with the `:child` field for all nodes in the subtree. 
"""
function istreeconsistent(node::T) where T <: AbstractTreeNode
    for n in getnodelist(node)
        for c in n.children
            if c.parent != n
                return false
            end
        end
    end
    return true
end

import Base.==
"""
    ==(t1, t2)
	Defaults to pointer equality
"""
function ==(t1::T, t2::T) where T <: AbstractTreeNode
	return ===(t1, t2)
end


export deepequals
"""
	deepequals(t1, t2)

Checks whether two trees are equal by recursively calling this on all fields, except `:parent`, in order to prevent cycles.
In order to ensure that the `:parent` field is not hiding something different on both trees, ensure that each is consistent first (see: `istreeconsistent`).
"""
function deepequals(t1::T, t2::T) where T <: AbstractTreeNode
    equal_or_both_undef = (x,y,f) -> (!isdefined(x,f) && !isdefined(y,f)) || (isdefined(x,f) && isdefined(y,f) && deepequals(getfield(x, f), getfield(y, f)))
    return mapreduce(f->equal_or_both_undef(t1, t2, f), (x,y) -> x && y, filter(x-> x != :parent, collect(fieldnames(GeneralFelNode))))
end


function deepequals(t1, t2)
	return t1 == t2
end


export addchild
function addchild(parent::T, child::T) where T <: AbstractTreeNode #This needs a bang!
  if child.parent === nothing
    push!(parent.children, child)
    #child.parent = Union{T,Nothing}(parent) #Unsure about this change.
    child.parent = parent
  else
    children = child.parent.children
    index = findfirst(children, child)
    deleteat!(children, index)
    child.parent = Union{T,Nothing}(parent)
    push!(parent.children, child)
  end
end


export safe_addchild
function safe_addchild(parent::T, child::T) where T <: AbstractTreeNode
  if child.parent === nothing
    push!(parent.children, child)
    #child.parent = Union{T,Nothing}(parent) #Unsure about this change.
    child.parent = parent
  else
    if insubtree(parent, child)
      throw(ArgumentError("Cannot move node to a subtree of itself!"))
    else
      children = child.parent.children
      index = findfirst(children, child)
      deleteat!(children, index)
      child.parent = Union{T,Nothing}(parent)
      push!(parent.children, child)
    end
  end
end


function mergenodes(n1::AbstractTreeNode, n2::AbstractTreeNode)
    newnode = eltype(n1)()
    if n1.parent === nothing && n2.parent === nothing
        addchild(newnode,n1)
        addchild(newnode,n2)
    else
        throw(ArgumentError("You are trying to merge two nodes that already have parents."))
    end
    return newnode
end

export insubtree
function insubtree(node::T, subtree::T) where T <: AbstractTreeNode
  if subtree == node
    return true
  else
    for child in subtree.children
      if insubtree(node, child)
        return true
      end
    end
  end
  return false
end

#=
function getnewickhelper(node::T) where T <: AbstractTreeNode
  if length(node.children) == 0
    return string(node.name,":", node.branchlength)
  else
    ret = join(AbstractString[getnewickhelper(child) for child in node.children], ",")
    return string("(",ret,")",node.name,":", node.branchlength)
  end
end

export getnewick
function getnewick(node::T) where T <: AbstractTreeNode
  return string(getnewickhelper(node),";")
end
=#


function treefromnewickhelper(newick::AbstractString)
  startindex = 1
  endindex = length(newick)
  lm = match(r"[\)][^\)]*$", newick)

  tag = newick
  childstrings = AbstractString[]
  if lm !== nothing
    tag = newick[lm.offset+1:end]
    childstring = newick[startindex+1:lm.offset-1]

    child = ""
    a = 0
    for i=1:length(childstring)
      if childstring[i] == '('
        a += 1
        child = string(child,childstring[i])
      elseif childstring[i] == ')'
        a -= 1
        child = string(child,childstring[i])
      elseif childstring[i] == ',' && a == 0
        push!(childstrings, child)
        child = ""
      else
        child = string(child,childstring[i])
      end
    end
    if child != ""
      push!(childstrings, child)
    end
  end
  spl = split(tag,":")
  name = ""
  branchlength = 0.0
  if length(spl) > 0
    name = strip(spl[1])
  end
  if length(spl) > 1
    branchlength = parse(Float64,spl[2])
  end

  return childstrings,name,branchlength
end

export gettreefromnewick
"""
  gettreefromnewick(newick)

Returns an AbstractTree of specified type from `newick` string.
"""

export gettreefromnewick 
function gettreefromnewick(str, T::DataType; tagged = false, disable_binarize = false)
    currnode = T()
    i = 1

    str = collect(str)

    function try_apply_char_arr(node, char_arr)
        if node.name == ""
            #println(char_arr)
            node.name = join(char_arr)
            #println("naming a node: ", join(char_arr))
        else
            #println("yo!")
        end
        empty!(char_arr)
    end
    char_arr = []
    tag_dict = Dict{T,Vector{String}}()
    while i <= length(str)
        c = str[i]
        
        if c == '{' && tagged
            init_loc = (i += 1)
            while i <= length(str)
                #println(i)
                if (str[i] == '}')
                    break
                end 
                i += 1
            end
            tag_dict[currnode] = string.(split(join(str[init_loc:i-1]), ","))
            i += 1
        elseif c == '(' 
            #println("making new node")
            try_apply_char_arr(currnode, char_arr)
            newnode = T() 
            addchild(currnode, newnode)
            currnode = newnode
            i += 1
        elseif c == ')' 
            try_apply_char_arr(currnode, char_arr)
            currnode = currnode.parent
            i += 1

        elseif c == ':' 
            try_apply_char_arr(currnode, char_arr)
            valid_nums = ['1','2','3','4','5','6','7','8','9','0','-','.', 'E', 'e']
            init_loc = (i += 1)
            while i <= length(str)
                #println(i)
                if !(str[i] in valid_nums)
                    break
                end 
                i += 1
            end 
            currnode.branchlength = parse(Float64, join(str[init_loc:i-1]))

        elseif c == ',' 
            #println("making new node")
            try_apply_char_arr(currnode, char_arr)
            newnode = T() 
            addchild(currnode.parent, newnode)
            currnode = newnode
            i += 1
        elseif c == ';' 
            try_apply_char_arr(currnode, char_arr)
            return (tagged ? (currnode, tag_dict) : currnode)
        else
            push!(char_arr, c)
            #println(char_arr)
            i += 1
        end 
    end

	binarize!(currnode) 
    
    return (tagged ? (currnode, tag_dict) : currnode)
end

export better_newick_import
function better_newick_import(str, T::DataType; tagged = false, disable_binarize = false)
	@warn "better_newick_import has been renamed and deprecated. Use gettreefromnewick instead."
	return gettreefromnewick(str, T; tagged=tagged, disable_binarize=disable_binarize)
end


export prettyprintstring
function prettyprintstring(node::T, spaces::Int=0) where T <: AbstractTreeNode
  ret = string(repeat("----",spaces),"+++++ ", node.name, "(", node.branchlength,")", "\n")
  for child in node.children
      ret = string(ret, prettyprintstring(child,spaces+1))
  end
  return ret
end

export getnodelist
function getnodelist(node::T, nodelist::Array{T,1}=T[]) where T <: AbstractTreeNode
  push!(nodelist,node)
  for childnode in node.children #Fixing this to avoid implementing iterate(::GeneralFelNode)
    getnodelist(childnode,nodelist)
  end
  return nodelist
end

function getpattern(data::Array{Float64,3}, node::T, col::Int, pattern::Array{Int8,1}=Int8[]) where T <: AbstractTreeNode
  if isleafnode(node)
    for b in data[node.seqindex,col,:]
      push!(pattern, Int8(b))
    end
  else
    for childnode in node.children
      getpattern(data,childnode,col,pattern)
    end
  end

  return pattern
end


export treedepth
function treedepth(node::T) where T <: AbstractTreeNode
    if isleafnode(node)
        return 1
    else
        return maximum([treedepth(node.children[i]) for i in 1:length(node.children)]) + 1
    end
end

export getnonleaflist
function getnonleaflist(node::T, nonleaflist::Array{T,1}=T[]) where T <: AbstractTreeNode
    if !isleafnode(node)
        push!(nonleaflist, node)
    end
    for childnode in node.children
        getnonleaflist(childnode, nonleaflist)
    end
    return nonleaflist
end

export getleaflist
function getleaflist(node::T, leaflist::Array{T,1}=T[]) where T <: AbstractTreeNode
    if isleafnode(node)
        push!(leaflist, node)
    end
    for childnode in node.children
        getleaflist(childnode, leaflist)
    end
    return leaflist
end

export getdistfromroot
function getdistfromroot(node::T) where T <: AbstractTreeNode
    if node.parent === nothing
        return 0
    else
        return getdistfromroot(node.parent) + node.branchlength
    end
end


function shuffletree(tree::T) where T <: AbstractTreeNode
    for node in getnodelist(tree)
        if length(node.children) != 0
            shuffle!(node.children)
        end
    end
    return tree
end

function ladderize(tree::T) where T <: AbstractTreeNode
	newtree = deepcopy(tree)
	ladderize!(newtree)
	return newtree
end

function ladderize!(tree::T) where T <: AbstractTreeNode
    for node in getnodelist(tree)
        if length(node.children) != 0
            sort!(node.children, lt= (x,y)->length(getnodelist(x)) < length(getnodelist(y)))
        end
    end
end


function getorder(tree::T) where T <: AbstractTreeNode
    return [node.seqindex for node in getleaflist(tree)]
end


function binarize(tree::T) where T <: AbstractTreeNode
    nodes = getnodelist(tree)
    counter = 0
    for n in nodes
        while length(n.children) > 2
            c1 = pop!(n.children)
            c2 = pop!(n.children)
            counter +=1
            push!(n.children, T(0.0, "binarized_$counter"))
            n.children[end].children = [c1,c2]
            n.children[end].parent = Union{T,Nothing}(n)
        end
    end
end

export siblings
"""
  siblings(node)

Returns a vector of siblings of node.
"""
function siblings(node::AbstractTreeNode)
    if isroot(node)
        return Array{typeof(newt)}([])
    else
        return [ch for ch in node.parent.children if ch != node]
    end
end

export sibling_inds
"""
  sibling_inds(node)

Returns logical indices of the siblings in the parent's child's vector.
"""
function sibling_inds(node::AbstractTreeNode)
    if isroot(node)
        return [false]
    else
        return Ref(node) .!= node.parent.children
    end
end

export name2node_dict
"""
  name2node_dict(root)

Returns a dictionary of leaf nodes, indexed by node.name. Can be used to associate sequences with leaf nodes.
"""
function name2node_dict(root::AbstractTreeNode)
    leaf_dict = Dict()
    for node in getleaflist(root)
        leaf_dict[node.name] = node
    end
    return leaf_dict
end

function reorient!(node::AbstractTreeNode, new_parent::AbstractTreeNode, new_branch_length)
    # Remove the new parent from the children list
    node = node
    #deleteat!(node.children, findfirst(node.children, new_parent))
    node.children = [n for n in node.children if n != new_parent]
    # If this node has a parent, reorient the parent.
    if node.parent !== nothing
        parent = node.parent
        reorient!(parent, node, node.branchlength)
        # Add the parent as a child
        push!(node.children, parent)
    end
    # Set then new parent as the parent, with the new branch length
    node.parent = new_parent
    node.branchlength = new_branch_length
end

function recursive_reroot!(child_node::AbstractTreeNode; dist_above_child=(child_node.branchlength/2))
    if dist_above_child > child_node.branchlength
        print("This isn't going to work")
    end

    # Remembering stuff
    dist_below_parent = child_node.branchlength - dist_above_child
    old_parent = child_node.parent

    new_root = typeof(child_node)()
    child_node.branchlength = dist_above_child
    reorient!(old_parent, child_node, dist_below_parent)
    new_root.children = [child_node, old_parent]
    child_node.parent = new_root
    old_parent.parent = new_root
    collapse_single_parents(new_root)
    return new_root
end

function recursive_reroot(child_node::AbstractTreeNode; dist_above_child=(child_node.branchlength/2))
	child_node = deepcopy(child_node)
    newtree = recursive_reroot!(child_node, dist_above_child=dist_above_child)
    return newtree
end

function collapse_single_parents(node::AbstractTreeNode)
    new_children = Array{typeof(node),1}()
    for child in node.children
        while length(child.children) == 1
            old_length = child.branchlength
            child = child.children[1]
            child.branchlength += old_length
        end
        collapse_single_parents(child)
        push!(new_children, child)
    end
    node.children = new_children
end

#For getting newick strings from tree objects.
function getnewickhelper(node::AbstractTreeNode)
  if length(node.children) == 0
    return string(node.name,":", node.branchlength)
  else
    ret = join(AbstractString[getnewickhelper(child) for child in node.children], ",")
    return string("(",ret,")",node.name,":", node.branchlength)
  end
end

export newick
function newick(node::AbstractTreeNode)
  return string(getnewickhelper(node),";")
end

function reroot!(node; dist_above_child=0)
    if dist_above_child > node.branchlength
        println("This isn't going to work")
    end 
    
    #begin by inserting a new fictitious node at the right location in the branch above the current child
    new_root = typeof(node)()
    #set up fictitious node
    node.parent.children[findfirst(x-> x==node, node.parent.children)] = new_root
    new_root.parent = node.parent
    #update child node
    new_root.children = [node]
    node.parent = new_root
    new_root.branchlength, node.branchlength = node.branchlength-dist_above_child, dist_above_child
    
    #build stack of nodes going up to the old root
    stack = [new_root]
    curr = new_root
    while !(isroot(curr))
        curr = curr.parent
        push!(stack, curr)
    end
        
    #special case: old root
    curr = pop!(stack)
    #bifurcation - the root should be removed
    if length(curr.children) == 2
        child_ind = findfirst(x->x != stack[end], curr.children)
        push!(stack[end].children, curr.children[child_ind])
        curr.children[child_ind].parent = stack[end]
        curr.children[child_ind].branchlength = curr.children[child_ind].branchlength + stack[end].branchlength
    else
        #let's just put the root back onto the stack, as it gets handled like any other node if we it is a polytomy
        push!(stack, curr)
    end
    
    while length(stack) > 1
        curr = pop!(stack)
        deleteat!(curr.children, findfirst(x->x==stack[end], curr.children))
        curr.parent = stack[end]
        curr.branchlength = curr.parent.branchlength
        push!(stack[end].children, curr)
    end
    
    new_root.branchlength, node.branchlength = node.branchlength-dist_above_child, dist_above_child
	new_root.parent = nothing   
 
    return new_root
end

function reroot(node; dist_above_child=0)
	return reroot!(deepcopy(node); dist_above_child=0)
end


function node2dist_traversal(node::AbstractTreeNode, distvec::Vector{Float64}, current_dist::Float64, up_pass::Bool, node_dic)
    if up_pass
        if !isroot(node)
            node2dist_traversal(node.parent, distvec, current_dist + node.branchlength, true, node_dic)
            for sib in siblings(node)
                node2dist_traversal(sib, distvec, current_dist + node.branchlength + sib.branchlength, false, node_dic)
            end
        end
    else
        if isleafnode(node)
            distvec[node_dic[node]] = current_dist
        else
            for ch in node.children
                node2dist_traversal(ch, distvec, current_dist + ch.branchlength, false, node_dic)
            end
        end
    end
end

"""
    tree2distances(root::AbstractTreeNode)

Returns a distance matrix for all pairs of leaf nodes, and a node-to-index dictionary.
Be aware that this dictionary will break when any of the node content (ie. anything on the tree) changes.
"""
function tree2distances(root::AbstractTreeNode)
    node_dic = Dict()
    dim = 0
    for (i,n) in enumerate(getleaflist(root))
        node_dic[n] = i
        dim = i
    end
    distmat = zeros(dim,dim)
    for n in getleaflist(root)
        distvec = zeros(dim)
        node2dist_traversal(n, distvec, 0.0, true, node_dic)
        distmat[:,node_dic[n]] = distvec
    end
    return distmat,node_dic
end

function root2tips_traversal(node::AbstractTreeNode, heightvec::Vector{Float64}, current_height::Float64, node_dic)
    if !isleafnode(node)
        for ch in node.children
            root2tips_traversal(ch, heightvec, current_height + ch.branchlength, node_dic)
        end
    else
        heightvec[node_dic[node]] = current_height
    end
end

"""
    root2tips(root::AbstractTreeNode)

Returns a vector of root-to-tip distances, and a node-to-index dictionary.
Be aware that this dictionary will break when any of the node content (ie. anything on the tree) changes.
"""
function root2tip_distances(root::AbstractTreeNode)
    node_dic = Dict()
    dim = 0
    for (i,n) in enumerate(getleaflist(root))
        node_dic[n] = i
        dim = i
    end
    heightvec = zeros(dim)
    root2tips_traversal(root,heightvec,0.0,node_dic)
    return heightvec,node_dic
end

"""
    tree2distances(root::AbstractTreeNode)

Returns a distance matrix for all pairs of leaf nodes, and a node-to-index dictionary.
Be aware that this dictionary will break when any of the node content (ie. anything on the tree) changes.
"""
function tree2shared_branch_lengths(root::AbstractTreeNode)
    node_dic = Dict()
    dim = 0
    for (i,n) in enumerate(getleaflist(root))
        node_dic[n] = i
        dim = i
    end
    distmat = zeros(dim,dim)
    for n in getleaflist(root)
        distvec = zeros(dim)
        node2dist_traversal(n, distvec, 0.0, true, node_dic)
        distmat[:,node_dic[n]] = distvec
    end
    heightvec = zeros(dim)
    root2tips_traversal(root,heightvec,0.0,node_dic)
    return abs.(heightvec .+ heightvec' .- distmat),node_dic
end
