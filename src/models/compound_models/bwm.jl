export BWMModel
"""
    mutable struct BWMModel{M} <: DiscreteStateModel where M <: DiscreteStateModel

# Fields
- `models::Vector{<:M}`: A vector of models.
- `weights::Vector{Float64}`: A vector of weights.

# Description
Branch-wise mixture model.
!!! note
    `forward!` and `backward!` are currently only defined for `M<:PMatrixModel`.
"""
mutable struct BWMModel{M} <: DiscreteStateModel where M <: DiscreteStateModel
    models::Vector{<:M}
    weights::Vector{Float64}
end

"""
    BWMModel{M}(models::Vector{<:M}) where M <: DiscreteStateModel

Convenience constructor where the weights are uniform.
"""
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