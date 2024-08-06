mutable struct BWMModel <: DiscreteStateModel
    models::Vector{<:DiscreteStateModel}
    weights::Vector{Float64}
end

function backward!(
    dest::DiscretePartition, 
    source::DiscretePartition, 
    model::BWMModel, 
    node::FelNode
)
    P = sum([getPmatrix(m,node) for m in model.models] .* (model.weights))
    mul!(dest.state, P, source.state)
    dest.scaling .= source.scaling
end

function forward!(
    dest::DiscretePartition, 
    source::DiscretePartition, 
    model::BWMModel, 
    node::FelNode
)
    P = sum([getPmatrix(m,node) for m in model.models] .* model.weights)
    dest.state .= (source.state'*P)'
    dest.scaling .= source.scaling
end
 
function eq_freq(model::BWMModel)
    sum([eq_freq(m) for m in model.models] .* model.weights)
end