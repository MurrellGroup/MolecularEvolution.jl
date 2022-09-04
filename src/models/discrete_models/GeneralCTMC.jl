mutable struct GeneralCTMC <: DiscreteStateModel
    Q::Array{Float64,2}
    r::Float64
    function GeneralCTMC(Q::Array{Float64,2})
        new(Q,1.0)
    end
end

function backward!(dest::DiscretePartition,
        source::DiscretePartition,
        model::GeneralCTMC,
        node::GeneralFelNode)
    P = exp(model.Q .* model.r .* node.branchlength)
    mul!(dest.state, P, source.state)
    dest.scaling .= source.scaling
end

function forward!(dest::DiscretePartition,
        source::DiscretePartition,
        model::GeneralCTMC,
        node::GeneralFelNode)
    P = exp(model.Q .* model.r .* node.branchlength)
    dest.state .= (source.state'*P)'
    dest.scaling .= source.scaling
end