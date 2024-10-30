export AlwaysUpModel
mutable struct AlwaysUpModel{T} <: BranchModel where {T <: BranchModel}
    model::T
end

function backward!(
    dest::Partition,
    source::Partition,
    model::AlwaysUpModel,
    node::FelNode,
)
    backward!(dest, source, model.model, node)
end

function forward!(
    dest::Partition,
    source::Partition,
    model::AlwaysUpModel,
    node::FelNode,
)
    backward!(dest, source, model.model, node)
end

function eq_freq(model::AlwaysUpModel)
    return eq_freq(model.model)
end