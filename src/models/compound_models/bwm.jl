export BWMModel
mutable struct BWMModel{M} <: DiscreteStateModel where M <: DiscreteStateModel
    models::Vector{<:M}
    weights::Vector{Float64}
end

function BWMModel{M}(models::Vector{<:M}) where M <: DiscreteStateModel
    BWMModel{M}(models, sum2one(ones(length(models))))
end

function backward!(
    dest::DiscretePartition, 
    source::DiscretePartition, 
    model::BWMModel{<:PMatrixModel}, 
    node::FelNode
)
    P = sum([getPmatrix(m,node) for m in model.models] .* (model.weights))
    mul!(dest.state, P, source.state)
    dest.scaling .= source.scaling
end

function forward!(
    dest::DiscretePartition, 
    source::DiscretePartition, 
    model::BWMModel{<:PMatrixModel}, 
    node::FelNode
)
    P = sum([getPmatrix(m,node) for m in model.models] .* model.weights)
    dest.state .= (source.state'*P)'
    dest.scaling .= source.scaling
end
 
function eq_freq(model::BWMModel)
    sum([eq_freq(m) for m in model.models] .* model.weights)
end