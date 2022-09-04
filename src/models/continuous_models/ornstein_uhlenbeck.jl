#Define a ZeroDriftOrnsteinUhlenbeck process, with minimal functionality for simulation
#(ie. branch_prop_down and eq freqs)
mutable struct ZeroDriftOrnsteinUhlenbeck <: ContinuousStateModel
    mean::Float64
    volatility::Float64
    reversion::Float64
    function ZeroDriftOrnsteinUhlenbeck(mean::Float64,volatility::Float64,reversion::Float64)
        new(mean,volatility,reversion)
    end
end

function forward!(dest::GaussianPartition,
        source::GaussianPartition,
        model::ZeroDriftOrnsteinUhlenbeck,
        node::FelNode)
        dest.mean =  (e^(-node.branchlength*model.reversion))*(source.mean - model.mean) + model.mean
        dest.var = ((1 - e^(-2*node.branchlength*model.reversion))*model.volatility^2)/(2*model.reversion)
        dest.norm_const = source.norm_const
end

#Need to define backward! to use this for inference

function eq_freq_from_template(model::ZeroDriftOrnsteinUhlenbeck,partition_template::GaussianPartition)
    out_partition = deepcopy(partition_template)
    out_partition.mean = model.mean
    out_partition.var = (model.volatility^2)/(2*model.reversion)
    out_partition.norm_const = 0.0
    return out_partition
end
