export CovarionModel
"""
# Constructors
```julia
CovarionModel(models::Vector{<:DiscreteStateModel}, inter_transition_rates::Matrix{Float64})
CovarionModel(models::Vector{<:DiscreteStateModel}, inter_transition_rate::Float64)
```
# Description
The covarion model.
"""
mutable struct CovarionModel <: DiscreteStateModel
    model::DiscreteStateModel
    function CovarionModel(models::Vector{<:DiscreteStateModel}, inter_transition_rates::Matrix{Float64})
        indiv_model_size = size(models[1].Q,1)
        num_models = length(models)
        total_size = indiv_model_size * num_models
        new_q = zeros(total_size, total_size)
        for i in 1:num_models
            select_range = ((i-1)*indiv_model_size)+1:(i*indiv_model_size)
            new_q[select_range, select_range] .= models[i].Q
        end
        #set transitions between same character state across different models
        for c in 1:indiv_model_size
            for i in 1:size(inter_transition_rates,1)
                for j in 1:size(inter_transition_rates,2)
                    new_q[((i-1)*indiv_model_size)+c, ((j-1)*indiv_model_size)+c] += inter_transition_rates[i,j]
                end
            end
        end
        for i in 1:size(new_q,1)
            new_q[i,i] = 0
            new_q[i,i] = -sum(new_q[i,:])
        end
        new(typeof(models[1])(new_q))
    end
    function CovarionModel(models::Vector{<:DiscreteStateModel}, inter_transition_rate::Float64)
        CovarionModel(models, 
            fill(inter_transition_rate, length(models), length(models)) .- 
            (Array{Float64,2}(I, length(models), length(models)) .* (inter_transition_rate * length(models)))
        )
    end
end
 
export CovarionPartition
#Parts is a vector containing the various components
#Weights is a num_components by num_sites matrix
"""
    CovarionPartition(states,sites,models,t)

A partition for the [`CovarionModel`](@ref).
"""
mutable struct CovarionPartition <: DiscretePartition
    state::Array{Float64,2}
    states::Int
    sites::Int
    scaling::Array{Float64,1}
    models::Int
    t::DataType
 
    function CovarionPartition(states,sites,models,t)
        new(zeros(states*models,sites),states*models,sites,models,t)
    end 
end
 
#"up" refers to propagating up the tree: ie time running backwards.
function backward!(
    dest::CovarionPartition,
    source::CovarionPartition,
    model::CovarionModel,
    node::FelNode
)
    backward!(dest, source, model.model, node)
end

#"down" refers to propagating up the tree: ie time running forwards.
function forward!(
    dest::CovarionPartition,
    source::CovarionPartition,
    model::CovarionModel,
    node::FelNode
)
    forward!(dest, source, model.model, node)
end
#One should either implement this funtion, or one that returns an entire "root partition".
function eq_freq(model::CovarionModel)
    return eq_freq(model.model)
end
 
function obs2partition!(dest::CovarionPartition,seq::String)
    new_part = dest.t(dest.sites)
    obs2partition!(new_part,seq)
    for i in 1:dest.models
        dest.state[((i-1)*new_part.states)+1:i*new_part.states,:] .= new_part.state
    end
end
 
function partition2obs(part::CovarionPartition)
    new_part = part.t(part.sites)
    for i in 1:part.models
        new_part.state .+= part.state[((i-1)*new_part.states)+1:i*new_part.states,:]
    end
    partition2obs(new_part)
end