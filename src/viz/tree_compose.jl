#Some of this code (especially the excellent highlighter plotting) was written by Alec Pankow.

#NOTE: I rely a lot on dictionaries with nodes as a key.
#This means that you cannot mutate the node contents and expect these to work.

function tree_lines(root::FelNode;
    dot_size_dict = Dict(), dot_size_default = 0.0,
    dot_color_dict = Dict(), dot_color_default = "black",
    line_color_dict = Dict(), line_color_default = "black",
    label_color_dict = Dict(), label_color_default = "black",
    nodelabel_dict = Dict())

    branch_Y_axis_counter = 0
    lines = []
    dots = []
    nodelabels = []

    dot_sizes = []
    line_colors = []
    dot_colors = []
    label_colors = []

    node_position_dict = Dict()

    function tree_lines_traversal(node::FelNode,lines,depth::Real)
        start_depth = depth
        depth += node.branchlength
        end_depth = depth
        #println(node.nodeindex," ",start_depth," ",end_depth," ",node.branchlength)
        Y_height = 0
        if isleafnode(node)
            Y_height = branch_Y_axis_counter
            branch_Y_axis_counter += 1
            push!(lines,[(start_depth,Y_height),(end_depth,Y_height)])
        else
            height_vec = [tree_lines_traversal(child,lines,end_depth) for child in node.children]
            Y_height = mean(height_vec)
            push!(lines,[(start_depth,Y_height),(end_depth,Y_height)])
            push!(lines,[(end_depth,maximum(height_vec)),(end_depth,minimum(height_vec))])

            push!(line_colors, get(line_color_dict,node,line_color_default))


        end
        push!(dots,(end_depth,Y_height))
        push!(nodelabels,get(nodelabel_dict,node,node.name))
        push!(dot_sizes, get(dot_size_dict,node,dot_size_default))
        push!(dot_colors, get(dot_color_dict,node,dot_color_default))
        push!(line_colors, get(line_color_dict,node,line_color_default))
        push!(label_colors, get(label_color_dict,node,label_color_default))

        node_position_dict[node] = (end_depth,Y_height)

        return Y_height
    end

    tree_lines_traversal(root,lines,0)

    maxdepth = maximum([max(i[1][1],i[2][1]) for i in lines])
    maxY = maximum([max(i[1][2],i[2][2]) for i in lines])

    scaled_lines = [[(i[1][1]/maxdepth,i[1][2]/maxY),(i[2][1]/maxdepth,i[2][2]/maxY)] for i in lines]
    scaled_dots = [(i[1]/maxdepth,i[2]/maxY) for i in dots]
    for k in keys(node_position_dict)
        node_position_dict[k] = (node_position_dict[k][1]/maxdepth,node_position_dict[k][2]/maxY)
    end
    return scaled_lines,scaled_dots,dot_sizes,nodelabels,dot_colors,line_colors,label_colors,node_position_dict
end

export simple_tree_draw
"""
img = simple_tree_draw(tree::FelNode; canvas_width = 15cm, canvas_height = 15cm, line_color = "black", line_width = 0.1mm)

A line drawing of a tree with very few options.

img = simple_tree_draw(tree)
img |> SVG("imgout.svg",10cm, 10cm)
OR
using Cairo
img |> PDF("imgout.pdf",10cm, 10cm)

"""
function simple_tree_draw(tree::FelNode; canvas_width = 15cm, canvas_height = 15cm, line_color = "black", line_width = 0.1mm)
    set_default_graphic_size(canvas_width, canvas_height)
    line_array,dot_array,dot_sizes,nodelabels,dot_colors,line_colors,label_colors,node_position_dict = tree_lines(tree)
    return compose(context(0.00, 0.00, 1.0, 1.0, units=UnitBox(-0.05,-0.05,1.1,1.1)),(context(),line(line_array), stroke(line_color), linewidth(line_width)))
end

function split_range(lb,stride,ub)
    pre_range = [lb:stride:ub;]
    if pre_range[end] <= ub
        push!(pre_range,ub)
    end
    return [pre_range[i]:pre_range[i+1] for i in 1:length(pre_range)-1]
end

export tree_draw
#Tree draw doesn't allow the vertical lines either side of a node to be colored differently - this isn't ideal and it should be re-done to allow this.
"""
    tree_draw(tree::FelNode;
        canvas_width = 15cm, canvas_height = 15cm,
        stretch_for_labels = 2.0, draw_labels = true,
        line_width = 0.1mm, font_size = 4pt,
        min_dot_size = 0.00, max_dot_size = 0.01,
        line_opacity = 1.0,
        dot_opacity = 1.0,
        name_opacity = 1.0,
        horizontal = true,
        dot_size_dict = Dict(), dot_size_default = 0.0,
        dot_color_dict = Dict(), dot_color_default = "black",
        line_color_dict = Dict(), line_color_default = "black",
        label_color_dict = Dict(), label_color_default = "black",
        nodelabel_dict = Dict(),compose_dict = Dict()
        )

Draws a tree with a number of self-explanatory options.
Dictionaries that map a node to a color/size are used to control per-node plotting options.
"compose_dict" must be a FelNode->function(x,y) dictionary that returns a compose() struct.

Example using compose_dict
str_tree = "(((((tax24:0.09731668728575642,(tax22:0.08792233964843627,tax18:0.9210388482867483):0.3200367900275155):0.6948314526087965,(tax13:1.9977212308725611,(tax15:0.4290074347886068,(tax17:0.32928401808187824,(tax12:0.3860215462534818,tax16:0.2197134841232339):0.1399122681886174):0.05744611946245004):1.4686085778061146):0.20724159879522402):0.4539334554156126,tax28:0.4885576926440158):0.002162260013924424,tax26:0.9451873777301325):3.8695419798779387,((tax29:0.10062813251515536,tax27:0.27653633028085006):0.04262434258357507,(tax25:0.009345653929737636,((tax23:0.015832941547076644,(tax20:0.5550597590956172,((tax8:0.6649025646927402,tax9:0.358506423199849):0.1439516404012261,tax11:0.01995439013213013):1.155181296134081):0.17930021667907567):0.10906638146207207,((((((tax6:0.013708993438720255,tax5:0.061144001556547097):0.1395453591567641,tax3:0.4713722705245479):0.07432598428904214,tax1:0.5993347898257291):1.0588025698844894,(tax10:0.13109032492533992,(tax4:0.8517302241963356,(tax2:0.8481963081549965,tax7:0.23754095940676642):0.2394313086297733):0.43596704123297675):0.08774657269409454):0.9345533723114966,(tax14:0.7089558245245173,tax19:0.444897137240675):0.08657675809803095):0.01632062723968511,tax21:0.029535281963725537):0.49502691718938285):0.25829576024240986):0.7339777396780424):4.148878039524972):0.0"
newt = better_newick_import(str_tree, FelNode{Float64});
newt = ladderize(newt);
compose_dict = Dict()
for n in getleaflist(newt)
    #Replace the rand(4) with the frequencies you actually want.
    compose_dict[n] = (x,y)->pie_chart(x,y,sum2one(rand(4)),size = 0.03)
end
tree_draw(newt,draw_labels = false,line_width = 0.5mm, compose_dict = compose_dict)


img = tree_draw(tree)
img |> SVG("imgout.svg",10cm, 10cm)
OR
using Cairo
img |> PDF("imgout.pdf",10cm, 10cm)
"""
function tree_draw(tree::FelNode;
        canvas_width = 15cm, canvas_height = 15cm,
        stretch_for_labels = 2.0, draw_labels = true,
        line_width = 0.1mm, font_size = 4pt,
        min_dot_size = 0.00, max_dot_size = 0.01,
        line_opacity = 1.0,
        dot_opacity = 1.0,
        name_opacity = 1.0,
        horizontal = true,
        dot_size_dict = Dict(), dot_size_default = 0.0,
        dot_color_dict = Dict(), dot_color_default = "black",
        line_color_dict = Dict(), line_color_default = "black",
        label_color_dict = Dict(), label_color_default = "black",
        nodelabel_dict = Dict(),
        compose_dict = Dict()
        )

    line_array,dot_array,dot_sizes,nodelabels,dot_colors,line_colors,label_colors,node_position_dict = tree_lines(tree,dot_size_dict = dot_size_dict, dot_size_default = dot_size_default,
        dot_color_dict = dot_color_dict, dot_color_default = dot_color_default,
        line_color_dict = line_color_dict, line_color_default = line_color_default,
        label_color_dict = label_color_dict, label_color_default = label_color_default,
        nodelabel_dict = nodelabel_dict)

    if maximum(dot_sizes) > 0
        dot_sizes ./= maximum(dot_sizes)
    end
    dot_sizes .*= max_dot_size
    dot_sizes = clamp.(dot_sizes,min_dot_size,max_dot_size)

    dotx = [i[1] for i in dot_array]
    doty = [i[2] for i in dot_array]
    maxX = stretch_for_labels
    maxY = 1.1

    vertical = !horizontal
    if vertical
        textx = dotx
        texty = doty .-(dot_sizes .+ mean(dot_sizes) .+ 0.005)
        rot = Rotation(pi/2,0.5,0.5)
        textrot = Rotation.([pi/2 for i in 1:length(dotx)], dotx, doty)
    else
        textx = dotx .+ dot_sizes .+ mean(dot_sizes) .+ 0.005
        texty = doty
        rot = Rotation(0,0.5,0.5)
        textrot = Rotation.([0 for i in 1:length(dotx)], dotx, doty)
    end

    line_ranges = split_range(1,50,length(line_array))
    dot_ranges = split_range(1,50,length(dotx))
    text_ranges = split_range(1,50,length(nodelabels))

    compose_keys = keys(compose_dict)
    custom = [compose_dict[k](node_position_dict[k][1],node_position_dict[k][2]) for k in compose_keys]

    set_default_graphic_size(canvas_width, canvas_height)
    img = compose(context(0.00, 0.00, 1.0, 1.0, units=UnitBox(-0.05,-0.05,maxX,maxY),rotation=rot),
        ifelse(length(compose_keys)>0,custom,compose(context())),
        ifelse(draw_labels,[compose(context(),Compose.text(textx[ran], texty[ran], nodelabels[ran], [hleft], [vcenter], textrot[ran]),fill(label_colors[ran]),stroke(label_colors[ran]),linewidth(0.1mm),fontsize(font_size),fillopacity(name_opacity),strokeopacity(name_opacity)) for ran in text_ranges],[compose(context())]),
        [compose(context(),circle(dotx[ran], doty[ran], dot_sizes[ran]), fill(dot_colors[ran]), fillopacity(dot_opacity)) for ran in dot_ranges],
        [compose(context(),line(line_array[ran]), stroke(line_colors[ran]), linewidth(line_width),strokeopacity(line_opacity)) for ran in line_ranges]
        )
    return img
end

#= Example
str_tree = readlines("Datasets/Wertheim2011/AIV_neuroaminidase.tre")[1];
newt = better_newick_import(str_tree, FelNode{Float64});
newt = ladderize(newt);
label_color_dict = Dict()
for n in getleaflist(newt)
    if occursin("mallard",n.name)
        label_color_dict[n] = "red"
    end
    if occursin("duck",n.name)
        label_color_dict[n] = "orange"
    end
end
tree_draw(newt,dot_color_default = "blue",draw_labels = true,label_color_dict = label_color_dict)
=#

export pie_chart
function pie_chart(locx,locy,freqs; colors = ["red","orange","green","blue", "purple"],size = 0.05,opacity = 1.0)
    if length(colors) < length(locx)
        error("Need a longer color vector.")
    end
    breaks = vcat([0],cumsum(freqs)*1.9999*π)
    return compose(context(), arc([locx], [locy], [size], breaks[1:end-1],breaks[2:end], [true]), fill(colors),fillopacity(opacity))
end


#http://algo.uni-konstanz.de/publications/bbs-dpt-05.pdf
mutable struct node_data
    descendents::Int
end

mutable struct node_loc
    pos::Tuple{Float64,Float64}
    w::Float64
    t::Float64
end

function descendent_count(node::AbstractTreeNode,node_data_dict)
    lv = 0
    if isleafnode(node)
        lv = 1
    else
        for ch in node.children
            lv = lv + descendent_count(ch,node_data_dict)
        end
    end
    node_data_dict[node] = lv
    return lv
end

function radial_layout(node::AbstractTreeNode,node_loc_dict,node_data_dict, total_descendents)
    v = node_loc_dict[node]
    if !isroot(node)
        u = node_loc_dict[node.parent]
        v.pos = u.pos .+ (node.branchlength .* (cos(v.t + v.w/2),sin(v.t + v.w/2)))
    end
    ng = v.t
    child_compose = []
    for ch in node.children
        cw = (node_data_dict[ch]/total_descendents)*2*pi
        ct = ng
        ng += cw
        node_loc_dict[ch] = node_loc((-99.0,-99.0),cw,ct)
        radial_layout(ch,node_loc_dict,node_data_dict,total_descendents)
    end
end

function radial_lines(root::FelNode,node_loc_dict)
    lines = Array{Array{Tuple{Float64,Float64},1},1}([])
    for node in getnodelist(root)
        if !isroot(node)
            push!(lines,[node_loc_dict[node.parent].pos,node_loc_dict[node].pos])
        end
    end
    return lines
end

function bounding_box(line_array)
    maxX = -Inf
    minX = Inf
    maxY = -Inf
    minY = Inf
    for l in line_array
        maxX = max(maxX,max(l[1][1],l[2][1]))
        maxY = max(maxY,max(l[1][2],l[2][2]))
        minX = min(minX,min(l[1][1],l[2][1]))
        minY = min(minY,min(l[1][2],l[2][2]))
    end
    return minX,minY,maxX,maxY
end

#need to test this one
#export simple_radial_tree_plot
"""
    simple_radial_tree_plot(root::FelNode; canvas_width = 10cm, line_color = "black", line_width = 0.1mm)

Draws a radial tree. No frills. No labels. Canvas height is automatically determined to avoid distorting the tree.

newt = better_newick_import("((A:1,B:1,C:1,D:1,E:1,F:1,G:1):1,(H:1,I:1):1);", FelNode{Float64});
simple_radial_tree_plot(newt,line_width = 0.5mm,root_angle = 7/10)
"""
function simple_radial_tree_plot(root::FelNode; canvas_width = 10cm, line_color = "black", line_width = 0.1mm, root_angle = 0.0)
    node_data_dict = Dict()
    total_descendents = descendent_count(root,node_data_dict)
    node_loc_dict = Dict()
    node_loc_dict[root] = node_loc((0.0,0.0),2*pi,root_angle)
    radial_layout(root,node_loc_dict,node_data_dict,total_descendents)
    line_array = radial_lines(root,node_loc_dict)
    #These all conspire to distort your tree. Need to handle this properly.
    minX,minY,maxX,maxY = bounding_box(line_array)
    maxX,maxY = maxX-minX,maxY-minY #These are like width and height
    #= If we want a square canvas
    minX = min(minX,minY)
    minY = minX
    maxX = max(maxX,maxY)
    maxY = maxX
    square_canvas_size = canvas_width
    set_default_graphic_size(square_canvas_size, square_canvas_size)
    =#
    set_default_graphic_size(canvas_width, canvas_width*(maxY/maxX))

    line_ranges = split_range(1,50,length(line_array))
    compose(context(0.00, 0.00, 1.0, 1.0, units=UnitBox(minX,minY,maxX,maxY)),
        [compose(context(),
                line(line_array[ran])
                , stroke(line_color),
                linewidth(line_width)
                ) for ran in line_ranges])
end


@require MultivariateStats="6f286f6a-111f-5878-ab1e-185364afe411" begin
    using .MultivariateStats

    export MDS_color_dict

    """
        MDS_color_dict(newt::AbstractTreeNode)

    Makes a node color dictionary which can be used to color tree leaves.
    d = MDS_color_dict(newt)
    tree_draw(newt,line_width = 0.5mm,label_color_dict = d)
    """
    function MDS_color_dict(newt::AbstractTreeNode)
        dis,dic = tree2distances(newt);
        #mds = classical_mds(dis,3);
        mds = transform(fit(MDS, dis, maxoutdim = 3, distances = true))
        PCA1 = mds[1,:]
        PCA2 = mds[2,:]
        PCA3 = mds[3,:]
        min1,max1 = minimum(PCA1),maximum(PCA1)
        min2,max2 = minimum(PCA2),maximum(PCA2)
        min3,max3 = minimum(PCA3),maximum(PCA3)
        redscale = linear_scale.(PCA1,min1,max1,0,1)
        bluescale = linear_scale.(PCA2,min2,max2,0,1)
        greenscale = linear_scale.(PCA3,min3,max3,0,1)
        color_function(i) = RGBA(greenscale[i],redscale[i],bluescale[i],1)
        colorscale = [color_function(i) for i in 1:length(redscale)];
        label_color_dict = Dict()
        leaflist = getleaflist(newt)
        for i in 1:length(leaflist)
            label_color_dict[leaflist[i]] = colorscale[i]
        end
        return label_color_dict
    end
    
end

export discrete_name_color_dict

"""
    discrete_name_color_dict(newt::AbstractTreeNode,tag_func; rainbow = false, scramble = false, darken = true, col_seed = nothing)

Takes a tree and a tag_func, which converts the leaf label into a category (ie. there should be <20 of these), and returns a color dictionary
that can be used to color the leaves or bubbles.

Example tag_func:
    function tag_func(nam::String)
        return split(nam,"_")[1]
    end

For prettier colors, but less discrimination: rainbow = true
To randomize the rainbow color assignment: scramble = true
col_seed is currently set to white, and excluded from the list of colors, to make them more visible.

Consider making your own version of this function to customize colors as you see fit.

Example use:
num_leaves = 50
Ne_func(t) = 1*(e^-t).+5.0
newt = sim_tree(num_leaves,Ne_func,1.0,nstart = rand(1:num_leaves));
newt = ladderize(newt)
tag_func(nam) = mod(sum(Int.(collect(nam))),7)
dic = discrete_name_color_dict(newt,tag_func,rainbow = true);
tree_draw(newt,line_width = 0.5mm,label_color_dict = dic)
"""
function discrete_name_color_dict(newt::AbstractTreeNode,tag_func; rainbow = false, scramble = false, darken = true, col_seed = RGB(1,1,1))
    tags = [tag_func(n.name) for n in getleaflist(newt)]

    uni_tags = union(tags)


    if rainbow
        if scramble
            uni_tags = sample(uni_tags,length(uni_tags),replace=false)
        end
        tag_dict = Dict()
        for i in 1:length(uni_tags)
            hue = (255((i-1)/(length(uni_tags))))
            tag_dict[uni_tags[i]] = HSVA(hue,1,1-(mod(i,2)*0.2),1)
            #print(hue, " ")
        end
    else
        cols = distinguishable_colors(length(uni_tags)+1, [col_seed])[2:end]
        tag_dict = Dict()
        for i in 1:length(uni_tags)
            tag_dict[uni_tags[i]] = cols[i]
        end
    end

    label_color_dict = Dict()
    leaflist = getleaflist(newt)
    for i in 1:length(leaflist)
        label_color_dict[leaflist[i]] = tag_dict[tag_func(leaflist[i].name)]
    end
    return label_color_dict
end

"""
    draw_example_tree(num_leaves = 50)

Draws a tree and shows the code that draws it.
"""
function draw_example_tree(;num_leaves = 50)
    Ne_func(t) = 1*(e^-t).+5.0
    newt = sim_tree(num_leaves,Ne_func,1.0,nstart = rand(1:num_leaves));
    newt = ladderize(newt)

    println(
    "
#Assuming \"newt\" is the root of the tree
d = MDS_color_dict(newt)
dot_size_dict = Dict()
for n in getleaflist(newt)
    dot_size_dict[n] = sqrt(rand())
end
tree_draw(newt,line_width = 0.5mm,label_color_dict = d,
    dot_size_dict = dot_size_dict, max_dot_size = 0.03,
    dot_color_dict = d, canvas_width = 10cm, canvas_height = 10cm)
    "
    )

    d = MDS_color_dict(newt)
    dot_size_dict = Dict()
    for n in getleaflist(newt)
        dot_size_dict[n] = sqrt(rand())
    end
    tree_draw(newt,line_width = 0.5mm,label_color_dict = d,dot_size_dict = dot_size_dict, max_dot_size = 0.03, dot_color_dict = d, canvas_width = 10cm, canvas_height = 10cm)
end

"""
    promote_internal(tree::FelNode)

Creates a new tree similar to the given tree, but with 'dummy' leaf nodes (w/ zero branchlength)
representing each internal node (for drawing / evenly spacing labels internal nodes).
"""
function promote_internal(tree::FelNode)
    newt = deepcopy(tree)
    for node in getnonleaflist(newt)
        ch = deepcopy(node)
        ch.children = []
        ch.parent = node
        ch.branchlength = 0
        if length(node.children) == 2
            # insert dummy in the middle
            insert!(node.children, 2, ch)
        else
            push!(node.children, ch)
        end
    end
#     MolecularEvolution.binarize(newt)
    return newt
end



"""
    highlight_seq_draw(x, y, str::AbstractString, region, basecolor, hicolor; fontsize=8pt, posx=hcenter, posy=vcenter)

Draw a sequence, highlighting the sites given in `region`.
This can be used along with `compose_dict` for drawing sequences at nodes in a tree (see `tree_draw`).
Returns a Compose container.
"""
function highlight_seq_draw(x, y, str::AbstractString, region, basecolor, hicolor; font_size=8pt, posx=hcenter, posy=vcenter)
    letters = []
    regs = Set(region)
    for i in 1:length(str)
        color = i in regs ? hicolor : basecolor
        push!(letters, (context(), text(i*font_size*0.75, 0, str[i], posx, posy), stroke(color)))
    end
    compose(context(x,y), letters..., fontsize(font_size))
end

const NUC_COLORS = Dict('A' => "red", 'G' => "gold", 'C' => "green3", 'T' => "blue")

const AA_COLORS = Dict(
    'D' => Compose.RGBA{Float64}(230,10,10,1), 'E' => Compose.RGBA{Float64}(230,10,10,1),
    'C' => Compose.RGBA{Float64}(230,230,0,1), 'M' => Compose.RGBA{Float64}(230,230,0,1),
    'R' => Compose.RGBA{Float64}(20,90,255,1), 'K' => Compose.RGBA{Float64}(20,90,255,1),
    'T' => Compose.RGBA{Float64}(250,150,0,1), 'S' => Compose.RGBA{Float64}(250,150,0,1),
    'F' => Compose.RGBA{Float64}(50,50,170,1), 'Y' => Compose.RGBA{Float64}(50,50,170,1),
    'N' => Compose.RGBA{Float64}(0,220,220,1), 'Q' => Compose.RGBA{Float64}(0,220,220,1),
    'G' => Compose.RGBA{Float64}(235,235,235,1),
    'L' => Compose.RGBA{Float64}(15,130,15,1), 'I' => Compose.RGBA{Float64}(15,130,15,1), 'V' => Compose.RGBA{Float64}(15,130,15,1),
    'A' => Compose.RGBA{Float64}(200,200,200,1),
    'W' => Compose.RGBA{Float64}(180,90,180,1),
    'H' => Compose.RGBA{Float64}(130,130,210,1),
    'P' => Compose.RGBA{Float64}(220,150,130,1),
)

"""
    colored_seq_draw(x, y, str::AbstractString; color_dict=Dict(), font_size=8pt, posx=hcenter, posy=vcenter)

Draw an arbitrary sequence.
`color_dict` gives a mapping from characters to colors (default black).
Default options for nucleotide colorings and amino acid colorings are given in the constants `NUC_COLORS` and `AA_COLORS`.
This can be used along with `compose_dict` for drawing sequences at nodes in a tree (see `tree_draw`).
Returns a Compose container.
"""
function colored_seq_draw(x, y, str::AbstractString; color_dict=Dict(), font_size=8pt, posx=hcenter, posy=vcenter)
    letters = []
    for i in 1:length(str)
        color = get(color_dict, str[i], "black")
        push!(letters, (context(), text(i*font_size*0.77, 0, str[i], posx, posy), stroke(color)))
    end
    compose(context(x,y), letters..., fontsize(font_size))
end

function dashed_line(x, y)
    ln = vcat([[(s,y),(s+0.02,y)] for s in x+0.01:0.03:1.06], [[(1.07,y),(1.09,y)]])
    return compose(context(), line(ln), linewidth(.05), stroke("black"))
end

function draw_partition(x, y, part::DiscretePartition; region=[], color_dict=Dict(), font_size=6pt)
    ns = partition2obs(part)
    if isempty(region)
        region = 1:length(ns)
    end
    return colored_seq_draw(1.1,y-0.05,ns[region],color_dict=color_dict,font_size=font_size,posy=vtop)
end

#color palettes for highlighter plots. I use an array of pairs
#to preserve color order

#cycle of AA colors, ordered by chemical properties
#added gap(-) and stop(*)
AA_HIGHLIGHT_COLORS = [
    'Q' => RGB(1.0,0.0,0.8),
    'E' => RGB(1.0,0.0,0.4),
    'D' => RGB(1.0,0.0,0.0),
    'S' => RGB(1.0,0.2,0.0),
    'T' => RGB(1.0,0.4,0.0),
    'G' => RGB(1.0,0.6,0.0),
    'P' => RGB(1.0,0.8,0.0),
    'C' => RGB(1.0,1.0,0.0),
    'A' => RGB(0.8,1.0,0.0),
    'V' => RGB(0.6,1.0,0.0),
    'I' => RGB(0.4,1.0,0.0),
    'L' => RGB(0.2,1.0,0.0),
    'M' => RGB(0.0,1.0,0.0),
    'F' => RGB(0.0,1.0,0.4),
    'Y' => RGB(0.0,1.0,0.8),
    'W' => RGB(0.0,0.8,1.0),
    'H' => RGB(0.0,0.4,1.0),
    'R' => RGB(0.0,0.0,1.0),
    'K' => RGB(0.4,0.0,1.0),
    'N' => RGB(0.8,0.0,1.0),
    '-' => "grey60",
    '*' => "black"
];

#NT colors with gap(-)
NT_HIGHLIGHT_COLORS = [
      'A' => "salmon",
      'G' => "gold",
      'T' => "lightgreen",
      'C' => "lightblue",
      '-' => "grey60",
];

"""
    get_max_depth(node,depth::Real)
Return the maximum depth of all children starting from the
indicated node.
"""
function get_max_depth(node,depth::Real)
    depth += node.branchlength
    if isleafnode(node)
        max_depth = depth
    else
        depth_vec = [get_max_depth(child,depth) for child in node.children]
        max_depth = maximum(depth_vec)
    end
    return max_depth
end

"""
    get_highlighter_legend(legend_colors)
Returns a Compose object given an input dictionary or pairs mapping
characters to colors.
"""
function get_highlighter_legend(legend_colors)
    width = 8
    spacing = (0.8/(maximum(mod.(0:(length(legend_colors) - 1),width))))
    coords = [(floor((i)/width; digits = 0)*0.6cm,spacing*(mod((i),width)) + 0.05) for i in 0:(length(legend_colors)-1)]
    x = [c[2] for c in coords]
    y = [c[1] for c in coords]
    compose(context(), fontsize(3),
        (context(), Compose.text(x .- 0.04, y .+ 0.3cm, ["$(k)" for (k,v) in legend_colors]), fill("black"), fillopacity(1.0)),
        (context(), rectangle(x,y,[0.4cm],[0.4cm]), fill([v for (k,v) in legend_colors])),
    )
end

#Return Compose object of the highlighted sequence for use in
#MolecularEvolution.tree_draw()
function get_highlighted_seq(x,y; lw = 0.5mm, node = node, seqnames = seqnames, seqs = seqs,
    tick_height = 1/(2.2*length(seqs)), master = master, highlighter_start = 1.1, highlighter_width = 1, color_dict = NT_HIGHLIGHT_COLORS)

    seq = seqs[findfirst(node.name .== seqnames)];
    mask = [c for c in uppercase(seq)] .!= [c for c in uppercase(master)]
    pos = collect(1:length(seq))[mask]
    chars = [c for c in seq][mask]
    color_array = [color_dict[c] for c in chars]
    point_array = [((ix/length(seq)*highlighter_width + highlighter_start,y-tick_height),(ix/length(seq)*highlighter_width + highlighter_start,y+tick_height)) for ix in pos]
    if length(point_array) > 0
        ranges = split_range(1,50,length(point_array))
    else
        ranges = [:]
    end
    #get composition
    return compose(context(0w,-0.05,1w,1.1),
        [compose(context(), line(point_array[r]), stroke(color_array[r]), linewidth(lw)) for r in ranges],
        (context(), line([((highlighter_start,y),(highlighter_width + highlighter_start,y))]), stroke("grey60"), linewidth(0.1mm)))
end

export highlighter_tree_draw

"""
    highlighter_tree_draw(tree, ali_seqs, seqnames, master;
        highlighter_start = 1.1, highlighter_width = 1,
        coord_width = highlighter_start + highlighter_width + 0.1,
        scale_length = nothing, major_breaks = 1000, minor_breaks = 500,
        tree_args = NamedTuple[], legend_padding = 0.5cm, legend_colors = NUC_colors)

Draws a combined tree and highlighter plot. The vector of seqnames must match
the node names in `tree`.

kwargs:
- tree_args: kwargs to pass to `tree_draw()`
- legend_colors: Mapping of characters to highlighter colors (default NT_colors)
- scale_length: Length of the scale bar
- highlighter_start: Canvas start for the highlighter panel
- highlighter_width: Canvas width for the highlighter panel
- coord_width: Total width of the canvas
- major_breaks: Numbered breaks for sequence axis
- minor_breaks: Ticks for sequence axis
"""
function highlighter_tree_draw(tree, ali_seqs, seqnames, master; highlighter_start = 1.1, highlighter_width = 1, coord_width = highlighter_start + highlighter_width + 0.1, scale_length = nothing,
    major_breaks = 1000, minor_breaks = 500, tree_args = NamedTuple[], legend_padding = 0.5cm, legend_colors = NT_HIGHLIGHT_COLORS, tick_width = 0.5mm)

    #parse taxa and construct highlighter
    compose_dict = Dict()
    for n in getleaflist(tree)
        compose_dict[n] = (x,y) -> get_highlighted_seq(x,y;node=n,highlighter_start=highlighter_start,highlighter_width=highlighter_width,master=master,
            seqnames=seqnames,seqs=ali_seqs,color_dict=Dict(legend_colors),lw=tick_width)
    end

    #legend
    legend = get_highlighter_legend(legend_colors)

    #axis
    seqlen = length(ali_seqs[1])
    start = (highlighter_start + 0.05)/coord_width #necessary to match `tree_draw()` start
    points_array = [((0,4.25mm),(1,4.25mm))]
    ticks_array = [((i/seqlen,4mm),(i/seqlen,6mm)) for i in 0:minor_breaks:seqlen]
    line_array = vcat(points_array, ticks_array)
    ax = compose(context(), line(line_array), linewidth(0.5mm),
    (context(), Compose.text([i/seqlen for i in 0:major_breaks:seqlen], [3mm], ["$(i)" for i in 0:major_breaks:seqlen], [Compose.hcenter]), linewidth(0.1mm)))

    #scale bar
    max_depth = get_max_depth(tree,0)
    if isnothing(scale_length) scale_length = round(max_depth/5,sigdigits=1) end
    scale_start = 0.5 - (scale_length/(2*max_depth))
    scale_end = 0.5 + (scale_length/(2*max_depth))
    scale_array = [((scale_start),1mm),((scale_end),(mm))]

    #draw
    img = tree_draw(tree; tree_args..., stretch_for_labels = coord_width, compose_dict = compose_dict)

    return compose(context(),
        (context(0,0.2cm,1,1h-legend_padding), img), #tree + highlighter
        (context(0.05/coord_width,1h-legend_padding,1/coord_width,1h), #scale bar
            line(scale_array), stroke("black"), linewidth(0.5mm),
            (context(), Compose.text(0.5,0,"$(scale_length)", Compose.hcenter), linewidth(0.1mm), fontsize(3), stroke("black"))),
        (context(start, 1h-legend_padding,1-start,1cm), legend), #legend
        (context(start,0,(highlighter_width)/coord_width,1h), ax, fontsize(3), stroke("black")) #axis
    )
end
