mutable struct GeneralCTMC <: DiscreteStateModel
    Q::Array{Float64,2}
    r::Float64
    function GeneralCTMC(Q::Array{Float64,2})
        new(Q, 1.0)
    end
end

getPmatrix(model::GeneralCTMC, node::FelNode) = exp(model.Q .* model.r .* node.branchlength)

"""
    backward!(dest::Partition, source::Partition, model::BranchModel, node::FelNode)

Propagate the source partition backwards along the branch to the destination partition, under the model.
Note: You should overload this for your own BranchModel types.
"""
function backward!(
    dest::DiscretePartition,
    source::DiscretePartition,
    model::GeneralCTMC,
    node::FelNode,
)
    P = getPmatrix(model, node)
    mul!(dest.state, P, source.state)
    dest.scaling .= source.scaling
end

"""
    forward!(dest::Partition, source::Partition, model::BranchModel, node::FelNode)

Propagate the source partition forwards along the branch to the destination partition, under the model.
Note: You should overload this for your own BranchModel types.
"""
function forward!(
    dest::DiscretePartition,
    source::DiscretePartition,
    model::GeneralCTMC,
    node::FelNode,
)
    P = getPmatrix(model, node)
    dest.state .= (source.state' * P)'
    dest.scaling .= source.scaling
end
