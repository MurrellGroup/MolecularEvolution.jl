
#Try replace these with randexp from Random
function rand_exp(lambda)
    rand(Exponential(1 / lambda))
end
function rand_exp()
    rand(Exponential(1.0))
end

function inhomo_poisson_next_sample(
    rate_func,
    start_time;
    min_delta = 0.001,
    rate_to_delta_scale = 10.0,
)
    #Slightly odd algorithm. Generate a "unit wait time" from an exponential with mean=1.
    #Then let the actual "amount" of rate accumulate until you've reached this unit wait.
    #Uses some heuristics for discretizing the function based on the current rate.
    #To do: Implement a piece-wise "exponential growth" version of this, and the user can 
    #pass in time,value pairs that will be linearly interpolated in the log domain.
    t = start_time
    current_rate = rate_func(t)
    delta = min(min_delta, 1 / (current_rate * rate_to_delta_scale))
    unit_wait = rand_exp()
    while unit_wait > delta
        current_rate = rate_func(t)
        unit_wait = unit_wait - current_rate * delta
        t += delta
        delta = min(min_delta, 1 / (current_rate * rate_to_delta_scale))
    end
    return t + unit_wait * delta
end

function sim_tree_function(
    add_limit::Int,
    Ne_func,
    sample_rate_func;
    nstart = 1,
    time = 0.0,
    mutation_rate = 1.0,
)
    nodes = Set((node = FelNode(0.0, "tax$(i)"), time = time) for i = 1:nstart)
    node_counter = length(nodes) + 1
    last_event = 0.0
    while length(nodes) > 1 || node_counter < add_limit
        flag = false
        ln = length(nodes)
        if ln > 1
            next_coal = inhomo_poisson_next_sample(
                t -> 1 / (2 * Ne_func(t) / (ln * (ln - 1))),
                time,
            ) #I should probs ask someone about this...
        else
            next_coal = Inf
        end
        if node_counter < add_limit
            next_add = inhomo_poisson_next_sample(sample_rate_func, time)
        else
            next_add = Inf
        end
        time = min(next_coal, next_add)
        if next_add < next_coal
            push!(nodes, (node = FelNode(0.0, "tax$(node_counter)"), time = time))
            last_event = time
            node_counter += 1
        elseif length(nodes) > 1
            left, right = sample2_without_replacement(nodes)
            merged = (node = mergenodes(left.node, right.node), time = time)
            left.node.branchlength = (time - left.time) * mutation_rate
            right.node.branchlength = (time - right.time) * mutation_rate
            delete!(nodes, left)
            delete!(nodes, right)
            push!(nodes, merged)
            last_event = time
            flag = true
        end
    end
    tree = pop!(nodes).node
    ladderize!(tree)
    return tree
end

function sample2_without_replacement(S)
    @assert length(S) > 1
    x = rand(S)
    y = rand(S)
    while isequal(y, x)
        y = rand(S)
    end
    return x, y
end

export sim_tree
"""
    sim_tree(add_limit::Int,Ne_func,sample_rate_func; nstart = 1, time = 0.0, mutation_rate = 1.0, T = Float64)

Simulates a tree of type FelNode{T}. Allows an effective population size function (Ne_func),
as well as a sample rate function (sample_rate_func), which can also just be constants.

Ne_func(t) = (sin(t/10)+1)*100.0 + 10.0
root = sim_tree(600,Ne_func,1.0)
simple_tree_draw(ladderize(root))
"""
function sim_tree(
    add_limit::Int,
    Ne,
    sample_rate;
    nstart = 1,
    time = 0.0,
    mutation_rate = 1.0,
)
    if Ne isa Number
        Ne_func = x -> Ne
    else
        Ne_func = Ne
    end
    if sample_rate isa Number
        sample_rate_func = x -> sample_rate
    else
        sample_rate_func = sample_rate
    end
    sim_tree_function(
        add_limit + 1,
        Ne_func,
        sample_rate_func;
        nstart = nstart,
        time = time,
        mutation_rate = mutation_rate,
    )
end

"""
    sim_tree(;n = 10)

Simulates tree with constant population size.

"""
function sim_tree(; n = 10, mutation_rate = 0.01)
    sim_tree(n, 500.0, 10.0, nstart = n, mutation_rate = mutation_rate)
end


"""
    standard_tree_sim(ntaxa)

Simulates a tree with logistic population growth, under a coalescent model.
"""
function standard_tree_sim(ntaxa)
    n(t) = (10*ntaxa)/(1+exp(t-10))
    return sim_tree(ntaxa,n,ntaxa/5, mutation_rate = 0.05)
end


"""
    ladder_tree_sim(ntaxa)

Simulates a ladder-like tree, using constant population size but heterochronous sampling, under a coalescent model.
"""
function ladder_tree_sim(ntaxa)
    n(t) = ntaxa/10
    return sim_tree(ntaxa,n,1.0, mutation_rate = 0.1/sqrt(ntaxa))
end