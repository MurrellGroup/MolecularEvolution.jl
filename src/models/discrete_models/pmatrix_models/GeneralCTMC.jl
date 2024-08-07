mutable struct GeneralCTMC <: PMatrixModel
    Q::Array{Float64,2}
    r::Float64
    function GeneralCTMC(Q::Array{Float64,2})
        new(Q, 1.0)
    end
end

getPmatrix(model::GeneralCTMC, node::FelNode) = exp(model.Q .* model.r .* node.branchlength)
