#BM: This needs to be updated after the change that avoids identity!, and put into some sort of documentation.

#A "model" is anything that specifies how the message at a node gets transformed as it moves from the lower node to the upper node. This can sometimes include specifying the actual branch prop itself, but is often just a Q matrix.
#Should contain all code related to specific models.

#When you define a new model, you also have to define some functions (some can be excluded if you only want some of the functionality):
#backward!(dest::NewPartitionType,source::NewPartitionType,model::NewModelType,node::FelNode)
#This must mutate the dest to have the message at the top of the branch.
#forward!(dest::NewPartitionType,source::NewPartitionType,model::NewModelType,node::FelNode)
#This must mutate the dest to have the message at the bottom of the branch.
#eq_freq(model::NewModelType)
#Returns the equilibrium frequency of the model. Typically for reversible models when getting root freqs.

#If you define a new partition type, that isn't a kind of DiscretePartition, then you need to define:
#combine!(dest::NewPartitionType,src::NewPartitionType)
#Combines src and dest and sets into dest
#identity!(dest::NewPartitionType)
#Sets dest to the identity value such that combine!(dest, other) == other
#site_LLs(dest::MyNewPartition)
#Gets a site-wise log-likelihood score

#For reference, anything <: DiscretePartition must have:
#state::Array{Float64,2}
#states::Int
#sites::Int
#scaling::Array{Float64, 1}
#With the data stored like states(rows),sites(columns).


#NOTE: add errors to all of these empty functions so that they explain why certain core functions don't work.
#---------------DEFINING A NEW PARTITION AND A NEW MODEL-----------------------------
#Copy this entire block. Fill in the functions you need for your purposes.
#PARTITIONS.
#Often you want to inherit DiscretePartition behavior because the combine, scaling, etc functions will just work, but if not:
mutable struct MyNewPartition <: Partition #Or another abstract partition type. Most common is: DiscretePartition{T}
    state::Array{Float64,2}    #Common
    states::Int                #but
    sites::Int                 #not 
    scaling::Array{Float64,1} #required
    function MyNewPartition(states::Int, sites::Int) #Implement a constructor.
        new(zeros(states, sites), states, sites)
    end
end
function combine!(dest::MyNewPartition, src::MyNewPartition)
    error("combine!() not yet implemented for $(typeof(dest)). Required for inference.")
    #dest.state .= f(dest.state,sr.state) #f must add the signal from sr onto dest.
end
#=
function identity!(dest::MyNewPartition)
    error("identity!() not yet implemented for $(typeof(dest)). Required for inference.")
    #fill!(dest.state,1.0) #This must set dest to the unit value, such that combine!(dest, other) == other
end
=#
function site_LLs(dest::MyNewPartition)
    error("site_LLs() not yet implemented for $(typeof(dest)). Required for inference.")
    #return f(dest) #f must produce an array of site-wise log-likelihoods. This can be of length 1 if you have only one site.
end

#MODELS.
mutable struct MyNewModel <: BranchModel #(or <: DiscreteStateModel or w/e)
    r::Float64 #Add fields. "r" is a common one, but not required.
    function MyNewModel(r::Float64) #Implement a constructor.
        new(r::Real)
    end
end
#"backward" refers to propagating up the tree: ie time running backwards.
function backward!(
    dest::MyNewPartition,
    source::MyNewPartition,
    model::MyNewModel,
    node::FelNode,
)
    error(
        "backward!() not yet implemented for $(typeof(model)) and $(typeof(source)). Required for any inference.",
    )
    #dest.state .= f(source.state) #You need to implement f(). Try and minimize allocation when this happens.
end
#"forward" refers to propagating up the tree: ie time running forwards.
function forward!(
    dest::MyNewPartition,
    source::MyNewPartition,
    model::MyNewModel,
    node::FelNode,
)
    error(
        "forward!() not yet implemented for $(typeof(model)) and $(typeof(source)). Required for simulation, and any internal state inference or branch length optimization.",
    )
    #dest.state .= f(source.state) #You need to implement f(). Try and minimize allocation when this happens.
end
#One should either implement this funtion, or one that returns an entire "root partition".
function eq_freq(model::MyNewModel)
    error(
        "eq_freq() not yet implemented for $(typeof(model)). Required for automatically setting root parent messages.",
    )
    #return Vector{T}
end
#---------------END: Every time you define a new model, implement these functions-----------------------------
