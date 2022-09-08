#To be a subtype of a DiscretePartition, a Partition must have ONLY state, states, sites, and scaling.
#Is there a way to enforce this?
#If there are extra fields, the copy_partition_to! will just ignore them, which isn't good.

#Overloading the copy_partition_to! to avoid allocations.
function copy_partition_to!(dest::T, src::T) where {T<:DiscretePartition}
    dest.state .= src.state
    dest.states = src.states
    dest.sites = src.sites
    dest.scaling .= src.scaling
end

#I should add a constructor that constructs a DiscretePartition from an existing array.
mutable struct CustomDiscretePartition <: DiscretePartition
    state::Array{Float64,2}
    states::Int
    sites::Int
    scaling::Array{Float64,1}
end

CustomDiscretePartition(states, sites) =
    CustomDiscretePartition(zeros(states, sites), states, sites, zeros(sites))

function CustomDiscretePartition(freq_vec::Vector{Float64}, sites::Int64) #Add this constructor to all partition types
    state_arr = zeros(length(freq_vec), sites)
    state_arr .= freq_vec
    return CustomDiscretePartition(state_arr, length(freq_vec), sites, zeros(sites))
end

mutable struct NucleotidePartition <: DiscretePartition
    state::Array{Float64,2}
    states::Int
    sites::Int
    scaling::Array{Float64,1}
end

NucleotidePartition(sites) = NucleotidePartition(zeros(4, sites), 4, sites, zeros(sites))

function NucleotidePartition(freq_vec::Vector{Float64}, sites::Int64)
    @assert length(freq_vec) == 4
    state_arr = zeros(4, sites)
    state_arr .= freq_vec
    return NucleotidePartition(state_arr, 4, sites, zeros(sites))
end

mutable struct GappyNucleotidePartition <: DiscretePartition
    state::Array{Float64,2}
    states::Int
    sites::Int
    scaling::Array{Float64,1}
end

GappyNucleotidePartition(sites) = GappyNucleotidePartition(zeros(5, sites), 5, sites, zeros(sites))

function GappyNucleotidePartition(freq_vec::Vector{Float64}, sites::Int64)
    @assert length(freq_vec) == 5
    state_arr = zeros(5, sites)
    state_arr .= freq_vec
    return GappyNucleotidePartition(state_arr, 5, sites, zeros(sites))
end

mutable struct AminoAcidPartition <: DiscretePartition
    state::Array{Float64,2}
    states::Int
    sites::Int
    scaling::Array{Float64,1}
end

AminoAcidPartition(sites) = AminoAcidPartition(zeros(20, sites), 20, sites, zeros(sites))

function AminoAcidPartition(freq_vec::Vector{Float64}, sites::Int64)
    @assert length(freq_vec) == 20
    state_arr = zeros(20, sites)
    state_arr .= freq_vec
    return AminoAcidPartition(state_arr, 20, sites, zeros(sites))
end

mutable struct GappyAminoAcidPartition <: DiscretePartition
    state::Array{Float64,2}
    states::Int
    sites::Int
    scaling::Array{Float64,1}
end

GappyAminoAcidPartition(sites) = GappyAminoAcidPartition(zeros(21, sites), 21, sites, zeros(sites))

function GappyAminoAcidPartition(freq_vec::Vector{Float64}, sites::Int64)
    @assert length(freq_vec) == 21
    state_arr = zeros(21, sites)
    state_arr .= freq_vec
    return GappyAminoAcidPartition(state_arr, 21, sites, zeros(sites))
end

function combine!(dest::DiscretePartition, src::DiscretePartition)
    # Update the state
    dest.state .*= src.state

    # Re-scaling to prevent underflow
    dest.scaling .+= src.scaling
    new_scaling = 1 ./ sum(dest.state, dims = 1)[:]
    scale_cols_by_vec!(dest.state, new_scaling)
    dest.scaling .+= -log.(new_scaling)
end

function identity!(dest::DiscretePartition)
    fill!(dest.state, 1.0)
    fill!(dest.scaling, 0.0)
end

#Replace each column with a 1-of-k vector that is a categorical draw from that row.
function sample_partition!(partition::DiscretePartition)
    for i = 1:partition.sites
        partition.state[:, i] .= one_hot_sample(partition.state[:, i])
    end
end


function max_partition!(part::DiscretePartition)
    for s = 1:part.sites
        m = argmax(part.state[:, s])
        part.state[:, s] .= 0.0
        part.state[m, s] = 1.0
    end
end

# Right after a combine, we should always have site_LLs == scaling - since we rescale after combine!
# However, if the last op we did was a prop the site_LLs here would be wrong
# TODO: robustness to last op being a prop
function site_LLs(dest::DiscretePartition)
    return dest.scaling
end

#Recovering existing behavior for other models. Need to replace the behavior of the current eq_freq code with this more general code.
function eq_freq_from_template(
    model::DiscreteStateModel,
    partition_template::DiscretePartition,
)
    out_partition = deepcopy(partition_template)
    eq_freq_vec = eq_freq(model) #This will dispatch correctly for existing models.
    out_partition.state .= repeat(eq_freq_vec, outer = [1, out_partition.sites])
    return out_partition
end

#Requires eq_freq() to be defined for each model.
function equilibrium_message(
    model_vec::Vector{<:DiscreteStateModel},
    message_template::Vector{<:DiscretePartition},
)
    out_mess = deepcopy(message_template)
    for part = 1:length(message_template)
        eq_freq_vec = eq_freq(model_vec[part])
        out_mess[part].state .= repeat(eq_freq_vec, outer = [1, out_mess[part].sites])
    end
    return out_mess
end
