#Diagonalized CTMC - currently fails if complex eigenvalues
#Need to specify a Q matrix
#Init triggers diagonalization
#Propagation requires matrix product
"""
# Constructors
```julia
DiagonalizedCTMC(Q::Array{Float64,2})
DiagonalizedCTMC(Qpre::Array{Float64,2}, pi::Vector{Float64})
```
# Description
Takes in a Q matrix (which can be multiplied onto row-wise by `pi`) and diagonalizes it.
When computing ``e^{Q t}`` (for different ``t``s), we now only need to exponentiate the eigenvalues,
and perform the two change-of-basis matrix multiplications.
!!! warning
    Construction fails if `Q` has complex eigenvalues.
"""
mutable struct DiagonalizedCTMC <: PMatrixModel
    Q::Array{Float64,2}
    D::Vector{Float64}
    V::Array{Float64,2}
    Vi::Array{Float64,2}
    r::Float64

    function DiagonalizedCTMC(Q::Array{Float64,2})
        D1, V1 = eigen(Q)
        Vi1 = pinv(V1)
        new(Q::Array{<:Real,2}, D1, V1, Vi1, 1.0)
    end

    function DiagonalizedCTMC(Qpre::Array{Float64,2}, pi::Vector{Float64})
        Q = Qpre .* pi'
        for i = 1:size(Q)[1]
            Q[i, i] = 0.0
            Q[i, i] = -sum(Q[i, :])
        end
        D1, V1 = eigen(Q)
        Vi1 = pinv(V1)
        new(Q, D1, V1, Vi1, 1.0)
    end
end

function eq_freq(model::DiagonalizedCTMC)
    pos = argmax(model.D)
    return (model.V[:, pos]) .* model.Vi[pos, :]
end

#=
function P_from_diagonalized_Q(model::DiagonalizedCTMC, node::FelNode)
    return clamp.(
        model.V * Diagonal(exp.(model.D .* model.r .* node.branchlength)) * model.Vi,
        0.0,
        Inf,
    )
end
=#

@deprecate P_from_diagonalized_Q(model::DiagonalizedCTMC, node::FelNode) getPmatrix(model::DiagonalizedCTMC, node::FelNode)
function getPmatrix(model::DiagonalizedCTMC, node::FelNode)
    return clamp.(
        model.V * Diagonal(exp.(model.D .* model.r .* node.branchlength)) * model.Vi,
        0.0,
        Inf,
    )
end
