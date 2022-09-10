mutable struct SWMModel <: BranchModel
    models::Vector{<:BranchModel}
    function SWMModel(models::Vector{<:BranchModel})
        new(models)
    end
    #Convenience constructor: takes a single model and an array of rates
    function SWMModel(model::M, rs::Vector{Float64}) where {M <: BranchModel}
        models = (vcat(model, [deepcopy(model) for x in rs[2:end]]))
        for (i,m) in enumerate(models)
            m.r = rs[i]
        end
        new(models)
    end
end

#Parts is a vector containing the various components
#Weights is a num_components by num_sites matrix
mutable struct SWMPartition{PType} <: Partition where {PType <: Partition}
    parts::Vector{PType}
    weights::Any
     sites::Int
    states::Int
    models::Int
    #constructor where sizes are sufficient to initialize parts
    function SWMPartition{PType}(parts::Vector{PType}) where {PType <: Partition}
        new{PType}(parts, nothing, parts[1].sites, parts[1].states, length(parts))
    end
    function SWMPartition{PType}(part::PType, n_parts::Int) where {PType <: Partition}
        SWMPartition{PType}([deepcopy(part) for part in 1:n_parts])
    end
end

function combine!(dest::SWMPartition{PType},src::SWMPartition{PType}) where {PType<:Partition}
    for i in 1:length(dest.parts)
        combine!(dest.parts[i], src.parts[i])
    end
end

function identity!(dest::SWMPartition)
    for part in dest.parts
        identity!(part)
    end
end

function site_LLs(dest::SWMPartition)
    error("site_LLs not implemented for SWMPartition")
end

#BM: Look into this:
function reweight(temp_message::SWMPartition, weights::Vector{Float64})
    tot = zero_scaling(temp_message.parts[1])
    for i in 1:length(temp_message.parts)
        tot .+= weights[i] .* temp_message.parts[i].scaling[((i-1)*length(tot))+1:(i * length(tot))]
    end
    return tot
end

#"backward" refers to propagating up the tree: ie time running backwards.
function backward!(dest::SWMPartition{PType},
        source::SWMPartition{PType},
        model::SWMModel,
        node::FelNode) where {PType<:Partition}
    for i in 1:length(dest.parts)
        backward!(dest.parts[i], source.parts[i], model.models[i], node)
    end
end
#"forward" refers to propagating up the tree: ie time running forwards.
function forward!(dest::SWMPartition{PType},
        source::SWMPartition{PType},
        model::SWMModel,
        node::FelNode) where {PType<:Partition}
    for i in 1:length(dest.parts)
        forward!(dest.parts[i], source.parts[i], model.models[i], node)
    end
end
#One should either implement this funtion, or one that returns an entire "root partition".
function eq_freq(model::SWMPartition)
    error("eq_freq() not yet implemented for $(typeof(model)). Required for automatically setting root parent messages.")
    #return Vector{T}
end

function eq_freq_from_template(model::SWMModel,partition_template::SWMPartition{PType}) where {PType<:Partition}
    out_partition = deepcopy(partition_template)
    for i in 1:length(out_partition.parts)
        out_partition.parts[i] = eq_freq_from_template(model.models[i], out_partition.parts[i])
    end
    return out_partition
end

function sample_partition!(partition::SWMPartition{PType}) where {PType <: Partition}
    for p in partition.parts
        sample_partition!(p)
    end
end

function string2partition!(dest::SWMPartition{PType},seq::String) where {PType <: DiscretePartition}
    for x in dest.parts
        string2partition!(x, seq)
    end
end

function SWM_sample_down!(tree; swm_index = 1)
    W = Weights(tree.message[swm_index].weights)
    model_selections = [sample(1:tree.message[swm_index].models, W) for i in 1:tree.message[swm_index].sites]
    for x in getleaflist(tree)
        for m in x.message[swm_index].parts
            for s in 1:tree.message[swm_index].sites
                m.state[:,s] = x.message[swm_index].parts[model_selections[s]].state[:,s]
            end
        end
    end
end

function total_LL(message::SWMPartition)
    return sum(reweight(message, message.weights))
end