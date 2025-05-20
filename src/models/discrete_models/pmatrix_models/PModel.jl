#This is literally just a single P matrix. Maybe some uses, but likely for testing speed bounds
"""
    PModel(P::Array{Float64,2})

A P matrix.
"""
mutable struct PModel <: PMatrixModel
    P::Array{Float64,2}
end

getPmatrix(model::PModel, node::FelNode) = model.P