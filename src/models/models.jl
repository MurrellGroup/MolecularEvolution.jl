include("continuous_models/continuous_models.jl")
include("discrete_models/discrete_models.jl")
include("compound_models/swm.jl")
include("lazy_models/lazy_partition.jl")


#BM: New way - avoids user having to define identity! function for new partitions - need to do a performance comparison vs the old way.

#Fallback method. This should be overloaded to avoid allocations wherever performance requires it
function copy_partition_to!(dest::T, src::T) where {T<:Partition}
    for f in fieldnames(T)
        setfield!(dest, f, deepcopy(getfield(src, f)))
    end
end

#Example overloading for GaussianPartition:
#=
function copy_partition_to!(dest::GaussianPartition,src::GaussianPartition)
    dest.mean = src.mean
    dest.var = src.var
    dest.norm_const = src.norm_const
end
=#

#Fallback method. This should be overloaded to avoid deepcopy wherever performance requires it
function copy_partition(src::T) where {T <: Partition}
    return deepcopy(src)
end

#Example overloading for GaussianPartition:
#=
function copy_partition(src::GaussianPartition)
    return GaussianPartition(src.mean, src.var, src.norm_const)
end
=#

#Return a new instance of T with the same shape as partition_template. Here, we don't care about the values stored within the partition.
#Fallback method. This should be overloaded with usage of undef to increase performance
function partition_from_template(partition_template::T) where {T <: Partition}
    return copy_partition(partition_template)
end

#Example overloading for DiscretePartition:
#=
function partition_from_template(partition_template::T) where {T <: DiscretePartition}
    states, sites = partition_template.states, partition_template.sites
    return T(Array{Float64, 2}(undef, states, sites), states, sites, Array{Float64, 1}(undef, sites))
end
=#

#Note: not enforcing a return type causes some unnecesarry conversions
copy_message(msg::Vector{<:Partition}) = [copy_partition(x) for x in msg]

#This is a function shared for all models - perhaps move this elsewhere
function combine!(dest::T, source_arr::Vector{<:T}, wipe::Bool) where {T<:Partition}
    #Init to be equal to 1, then multiply everything on.
    if wipe
        copy_partition_to!(dest, source_arr[1])
        for i = 2:length(source_arr)
            combine!(dest, source_arr[i])
        end
    else
        for src in source_arr
            combine!(dest, src)
        end
    end
end
