#InterpolatedDiscreteModel works by storing a number of P matrices, and the "t" values to which they correspond
#For a requested t, the returned P matrix is interpolated between it's two neighbours

function check_eq_P(P)
    return maximum(std(P,dims = 1))
end

mutable struct InterpolatedDiscreteModel <: DiscreteStateModel
    tvec::Vector{Float64}
    Pvec::Array{Float64,3} #Now a tensor. Pvec[:,:,i] is the ith P matrix.
    #Generator must be something that returns a P matrix with T as the only argument
    function InterpolatedDiscreteModel(siz::Int64, generator, tvec::Vector{Float64})
        @assert tvec[1] == 0.0
        @assert issorted(tvec)
        Pvec = zeros(siz,siz,length(tvec))
        for (i,t) in enumerate(tvec)
            Pvec[:,:,i] .= generator(t)
        end
        cp = check_eq_P(Pvec[:,:,end])
        if cp > 10^-9
            @warn "Max std dev of last P matrix is $(cp). Far from equilibrium - extend range"
        end
        new(tvec,Pvec)
    end
    #Alternative constructor where you directly specify the tensor of P matrices for each "t"
    function InterpolatedDiscreteModel(Pvec::Array{Float64,3}, tvec::Vector{Float64})
        @assert tvec[1] == 0.0
        @assert issorted(tvec)
        cp = check_eq_P(Pvec[:,:,end])
        if cp > 10^-9
            @warn "Max std dev of last P matrix is $(cp). Far from equilibrium - extend range"
        end
        new(tvec,Pvec)
    end
end


#The keys in d must be the boundaries of the ranges for which we would index into
function range_index(v::Float64,ts::Vector{Float64})
    p = searchsortedlast(ts,v) # index of the last key less than or equal to v
    return p,p+1
end

function interp_weight(tup,p)
    @assert tup[1] <= p && p <= tup[2]
    w = (p - tup[1]) ./ (tup[2] - tup[1])
    return 1-w,w
end

#We could maybe hang a destintion "matrix" on the model, and this would store the
#interpolated P to that, saving new allocations.
function matrix_interpolate(t_query,ts,Ps)
    inds = range_index(t_query,ts)
    if inds[2] > length(ts)
        return Ps[:,:,end]
    end
    w = interp_weight((ts[inds[1]],ts[inds[2]]),t_query)
    approxP = w[1].*Ps[:,:,inds[1]] .+ w[2].*Ps[:,:,inds[2]]
    return approxP
end

function backward!(
        dest::DiscretePartition,
        source::DiscretePartition,
        model::InterpolatedDiscreteModel,
        node::FelNode)
    P = matrix_interpolate(node.branchlength, model.tvec, model.Pvec)
    mul!(dest.state, P, source.state)
    dest.scaling .= source.scaling
end

function forward!(
        dest::DiscretePartition, 
        source::DiscretePartition, 
        model::InterpolatedDiscreteModel, 
        node::FelNode)
    P = matrix_interpolate(node.branchlength, model.tvec, model.Pvec)
    dest.state .= (source.state'*P)'
    dest.scaling .= source.scaling
end

function eq_freq(model::InterpolatedDiscreteModel)
    model.Pvec[1,:,end]
end

#step: Higher numbers mean smaller jumps
#cap: After this many points it starts doubling
function t_sequence(t::Float64,n::Int64; step = 2 ,cap = n - 10) 
    ts = zeros(n)
    ts[1] = 0.0 #Note setting the first one to the zero
    ts[2] = t
    c = 2
    for i in 3:n
        ts[i] = ts[i-1]+ts[c]
        if mod(i,step)==0
            c += 1
        elseif i > cap
            c = i
        end
    end
    return ts
end

function matrix_sequence(Q::Array{Float64,2},t::Float64,n::Int64; step = 2 ,cap = n - 10) 
    P = exp(Q .* t)
    Ps = zeros(size(P)[1],size(P)[2],n) #Big stack of matrices
    Ps[:,:,1] .= Diagonal(ones(size(P)[1])) #Note setting the first one to the identity
    Ps[:,:,2] .= P
    c = 2
    for i in 3:n
        Ps[:,:,i] .= Ps[:,:,i-1]*Ps[:,:,c]
        if mod(i,step)==0
            c += 1
        elseif i > cap
            c = i
        end
    end
    return Ps
end

#This will take an existing InterpolatedDiscreteModel, and effectively scale all the "t" values in e^Qt.
function rescale!(m::InterpolatedDiscreteModel, factor::Float64)
    m.tvec .= m.tvec ./ factor
end
