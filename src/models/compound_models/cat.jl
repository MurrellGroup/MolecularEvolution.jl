mutable struct CATModel <: BranchModel
    models::Vector{<:BranchModel}
    function CATModel(models::Vector{<:BranchModel})
        new(models)
    end
end
 
mutable struct CATPartition{PType} <: Partition where {PType <: DiscretePartition}
    part_inds::Vector{Vector{Int}}
    parts::Vector{PType}
    #constructor where sizes are sufficient to initialize parts
    function CATPartition{PType}(part_inds::Vector{Vector{Int}}) where {PType <: Partition}
        new{PType}(part_inds, Vector{PType}([PType(length(x)) for x in part_inds]))
    end
    #constructor where someone passes in instantiated partitions as well
    function CATPartition{PType}(part_inds::Vector{Vector{Int}}, parts::Vector{PType}) where {PType <: Partition}
        new{PType}(part_inds, parts)
    end
end
 
function combine!(dest::CATPartition{PType},src::CATPartition{PType}) where {PType<:DiscretePartition}
    for i in 1:length(dest.parts)
        combine!(dest.parts[i], src.parts[i])
    end
end
 
function site_LLs(part::CATPartition{PType}) where {PType<:DiscretePartition}
    ret_val = zeros(sum(length.(p.part_inds)))
    for (i, inds) in enumerate(part.part_inds)
        ret_val[inds] = site_LLs(parts[i])
    end
end
 
function identity!(dest::CATPartition{PType}) where {PType<:DiscretePartition}
    for p in parts
        identity!(p)
    end
end
 
#"up" refers to propagating up the tree: ie time running backwards.
function backward!(
    dest::CATPartition{PType},
    source::CATPartition{PType},
    model::CATModel,
    node::FelNode
) where {PType<:DiscretePartition}
    for i in 1:length(dest.parts)
        backward!(dest.parts[i], source.parts[i], model.models[i], node)
    end
end
 
#"down" refers to propagating up the tree: ie time running forwards.
function forward!(
    dest::CATPartition{PType},
    source::CATPartition{PType},
    model::CATModel,
    node::FelNode
) where {PType<:DiscretePartition}
    for i in 1:length(dest.parts)
        forward!(dest.parts[i], source.parts[i], model.models[i], node)
    end
end
 
#One should either implement this funtion, or one that returns an entire "root partition".
function eq_freq(model::CATModel)
    error("eq_freq() not yet implemented for $(typeof(model)). Required for automatically setting root parent messages.")
    #return Vector{T}
end
 
function eq_freq_from_template(model::CATModel,partition_template::CATPartition{PType}) where {PType <: DiscretePartition}
    out_partition = deepcopy(partition_template)
    for i in 1:length(out_partition.parts)
        out_partition.parts[i] = eq_freq_from_template(model.models[i], out_partition.parts[i])
    end
    return out_partition
end
 
function sample_partition!(partition::CATPartition{PType}) where {PType <: DiscretePartition}
    for p in partition.parts
        sample_partition!(p)
    end
end
 
#BM: This will fail for eg. codon models, where 3 chars becomes a codon.
function obs2partition!(dest::CATPartition{PType},seq::String) where {PType <: DiscretePartition}
    for i in 1:length(dest.part_inds)
        obs2partition!(dest.parts[i], seq[dest.part_inds[i]])
    end
end
 
function partition2obs(part::CATPartition{PType}) where {PType <: DiscretePartition}
    out = Vector{Char}(undef, sum(length.(part.part_inds)))
    for i in 1:length(part.parts)
        out[part.part_inds[i]] .= collect(partition2obs(part.parts[i]))
    end
    return join(out)
end
 