include("models_update.jl")
include("AbstractUpdate.jl")
include("felsenstein.jl")
include("branchlength_optim.jl")
include("lls.jl")
include("nni_optim.jl")
include("root_optim.jl")
include("ancestors.jl")
include("generative.jl")

#Maybe we should use safepop! for LazyPartition too?
function safepop!(temp_messages::Vector{Vector{T}}, temp_message::Vector{T}) where T <: Partition
    return isempty(temp_messages) ? copy_message(temp_message) : pop!(temp_messages)
end