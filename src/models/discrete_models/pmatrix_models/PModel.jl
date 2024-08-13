#This is literally just a single P matrix. Maybe some uses, but likely for testing speed bounds
mutable struct PModel <: DiscreteStateModel
    P::Array{Float64,2}
end

getPmatrix(model::PModel) = model.P