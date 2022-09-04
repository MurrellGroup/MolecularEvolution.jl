#Diagonalized CTMC - currently fails if complex eigenvalues
#Need to specify a Q matrix
#Init triggers diagonalization
#Propagation requires matrix product
mutable struct DiagonalizedCTMC <: DiscreteStateModel
    Q::Array{Float64,2}
    D::Vector{Float64}
    V::Array{Float64,2}
    Vi::Array{Float64,2}
    r::Float64

    function DiagonalizedCTMC(Q::Array{Float64,2})
        D1,V1 = eigen(Q) 
        Vi1 = pinv(V1)
        new(Q::Array{<:Real,2},D1,V1,Vi1,1.0)
    end

    function DiagonalizedCTMC(Qpre::Array{Float64,2},pi::Vector{Float64})
        Q = Qpre .* pi'
        for i in 1:size(Q)[1]
            Q[i,i] = 0.0
            Q[i,i] = -sum(Q[i,:])
        end
        D1,V1 = eigen(Q) 
        Vi1 = pinv(V1)
        new(Q,D1,V1,Vi1,1.0)
    end
end

function eq_freq(model::DiagonalizedCTMC)
    pos = argmax(model.D)
    return (model.V[:,pos]).*model.Vi[pos,:]
end

function P_from_diagonalized_Q(model::DiagonalizedCTMC,node::GeneralFelNode)
    return clamp.(model.V*Diagonal(exp.(model.D.*model.r.*node.branchlength))*model.Vi,0.0,Inf)
end


#"backward" refers to propagating up the tree: ie time running backwards.
function backward!(dest::DiscretePartition,
        source::DiscretePartition,
        model::DiagonalizedCTMC,
        node::GeneralFelNode)

    #P = model.V*Diagonal(exp.(model.D.*model.r.*node.branchlength))*model.Vi
    P = P_from_diagonalized_Q(model,node)
    mul!(dest.state, P, source.state)
    dest.scaling .= source.scaling
end



#Model list should be a list of P matrices.
function forward!(dest::DiscretePartition,
        source::DiscretePartition,
        model::DiagonalizedCTMC,
        node::GeneralFelNode)

    P = P_from_diagonalized_Q(model,node)
    dest.state .= (source.state'*P)' #Perf check here?
    dest.scaling .= source.scaling
end