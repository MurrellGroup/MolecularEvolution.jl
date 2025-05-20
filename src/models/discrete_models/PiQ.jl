#F81 but for general dimensions.

#Propagation in O(States)
"""
# Constructors
```julia
PiQ(r::Float64,pi::Vector{Float64}; normalize=false)
PiQ(pi::Vector{Float64}; normalize=false)
```
# Description
The F81 substitution model, but for general dimensions.
https://www.diva-portal.org/smash/get/diva2:1878793/FULLTEXT01.pdf
"""
mutable struct PiQ <: DiscreteStateModel
    r::Float64
    pi::Vector{Float64}
    beta::Float64

    function PiQ(r::Float64,pi::Vector{Float64}; normalize=false)
        piNormed = sum2one(pi)
        beta = normalize ? 1/(1-sum(abs2.(piNormed))) : 1.0
        new(r,piNormed,beta)
    end
    function PiQ(pi::Vector{Float64}; normalize=false)
        PiQ(1.0,pi;normalize=normalize)
    end
end

function backward!(dest::DiscretePartition,
        source::DiscretePartition,
        model::PiQ,
        node::FelNode)
    pow = exp.(-model.beta*model.r*node.branchlength)
    c1 = ((1 - pow).*model.pi)
    c2 = (pow .+ c1)
    vsum = sum(source.state .* c1, dims=1)
    dest.state .= pow .* source.state .+ vsum
    dest.scaling .= source.scaling
end

function forward!(dest::DiscretePartition,
        source::DiscretePartition,
        model::PiQ,
        node::FelNode)
    #ToDo: check this is the same as F81 using the full matrix exponential.
    scals = sum(source.state,dims = 1)[:] #Requires V1.0 fix.
    pow = exp.(-model.beta*model.r*node.branchlength)
    c1 = ((1 - pow).*model.pi)
    c2 = (pow .+ ((1 - pow).*model.pi))
    dest.state .= (scals' .- source.state).*c1 .+ source.state.*c2
    dest.scaling .= source.scaling
end

function eq_freq(model::PiQ)
    return model.pi
end

export PiQ
