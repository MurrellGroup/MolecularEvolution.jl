#Model behavior
mutable struct BrownianMotion <: ContinuousStateModel
    mean_drift::Float64
    var_drift::Float64
    function BrownianMotion(mean_drift::Float64, var_drift::Float64)
        new(mean_drift, var_drift)
    end
end

function backward!(
    dest::GaussianPartition,
    source::GaussianPartition,
    model::BrownianMotion,
    node::FelNode,
)
    dest.mean = source.mean - node.branchlength * model.mean_drift
    dest.var = source.var + node.branchlength * model.var_drift
    dest.norm_const = source.norm_const
end

function forward!(
    dest::GaussianPartition,
    source::GaussianPartition,
    model::BrownianMotion,
    node::FelNode,
)
    dest.mean = source.mean + node.branchlength * model.mean_drift
    dest.var = source.var + node.branchlength * model.var_drift
    dest.norm_const = source.norm_const
end

#If you want to use a root prior, you can set these values explicitly.
function eq_freq_from_template(model::BrownianMotion, partition_template::GaussianPartition)
    out_partition = copy_partition(partition_template)
    out_partition.mean = 0.0
    out_partition.var = Inf
    out_partition.norm_const = 0.0
    return out_partition
end
