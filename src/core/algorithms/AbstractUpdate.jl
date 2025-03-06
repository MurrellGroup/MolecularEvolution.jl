"""
# Summary
`abstract type AbstractUpdate <: Function`

A callable type that typically takes `(tree::FelNode, models; partition_list=1:length(tree.message))`, updates `tree` and `models`, and returns the updated `tree` and `models`.
# Example
Define a new subtype, where `foo` and `bar` are arbitrary updating functions

```julia
struct MyUpdate <: AbstractUpdate end

function (update::MyUpdate)(tree::FelNode, models; partition_list=1:length(tree.message))
    tree, models = foo(tree, models, partition_list=partition_list)
    tree, models = BayesUpdate(nni=0)(tree, models, partition_list=partition_list)
    tree, models = bar(tree, models, partition_list=partition_list)
    return tree, models
end
```
See also: [`StandardUpdate`](@ref)
"""
abstract type AbstractUpdate <: Function end

"""
# Summary
`struct StandardUpdate <: AbstractUpdate`

A standard update can be a family of calls to [`nni_update!`](@ref), [`branchlength_update!`](@ref), [`root_update!`](@ref), and model updates.
# Constructor
    StandardUpdate(
        nni::Int,
        branchlength::Int,
        root::Int,
        models::Int,
        refresh::Bool,
        nni_selection::Function,
        branchlength_modifier::UnivariateModifier,
        root_update::RootUpdate,
        models_update::ModelsUpdate
    )

# Arguments
- `nni::Int`: the number of times to update the tree by `nni_update!`
- `branchlength::Int`: the number of times to update the tree by `branchlength_update!`
- `root::Int`: the number of times to update the tree by `root_update!`
- `models::Int`: the number of times to update the model
- `refresh::Bool`: whether to refresh the messages in tree between update operations to ensure message consistency
- `nni_selection::Function`: the function that selects between nni configurations
- `branchlength_modifier::UnivariateModifier`: the modifier to update a branchlength by `branchlength_update!`
- `root_update::RootUpdate`: updates the root by `root_update!`
- `models_update::ModelsUpdate`: updates the model parameters

See also: [`BayesUpdate`](@ref), [`MaxLikUpdate`](@ref)
"""
struct StandardUpdate <: AbstractUpdate
    nni::Int
    branchlength::Int
    root::Int
    models::Int
    refresh::Bool
    nni_selection::Function
    branchlength_modifier::UnivariateModifier
    root_update::RootUpdate
    models_update::ModelsUpdate
end

"""
    BayesUpdate(;
        nni = 1,
        branchlength = 1,
        root = 0,
        models = 0,
        refresh = false,
        branchlength_sampler::UnivariateSampler = BranchlengthSampler(
            Normal(0, 2),
            Normal(-1, 1),
        ),
        root_sampler::RootSample = StandardRootSample(1.0, 1),
        models_sampler::ModelsUpdate = StandardModelsUpdate()
    )

Convenience constructor for [`StandardUpdate`](@ref). The `nni_selection` is fixed to `softmax_sampler`. 
This constructor provides Bayesian updates by sampling from the posterior distribution.
"""
BayesUpdate(;
    nni::Int = 1,
    branchlength::Int = 1,
    root::Int = 0,
    models::Int = 0,
    refresh::Bool = false,
    branchlength_sampler::UnivariateSampler = BranchlengthSampler(
        Normal(0, 2),
        Normal(-1, 1),
    ),
    root_sampler = StandardRootSample(1e-2, 1),
    models_sampler::ModelsUpdate = StandardModelsUpdate(),
) = StandardUpdate(
    nni,
    branchlength,
    root,
    models,
    refresh,
    softmax_sampler,
    branchlength_sampler,
    root_sampler,
    models_sampler,
)

"""
    MaxLikUpdate(;
        nni = 1,
        branchlength = 1,
        root = 0,
        models = 0,
        refresh = false,
        branchlength_optimizer::UnivariateOpt = GoldenSectionOpt(),
        root_optimizer = StandardRootOpt(10),
        models_optimizer::ModelUpdate = StandardModelUpdate()
    )

Convenience constructor for [`StandardUpdate`](@ref). The `nni_selection` is fixed to `argmax`.
This constructor provides Maximum Likelihood updates by optimizing parameters.
"""
MaxLikUpdate(;
    nni::Int = 1,
    branchlength::Int = 1,
    root::Int = 0,
    models::Int = 0,
    refresh::Bool = false,
    branchlength_optimizer::UnivariateOpt = GoldenSectionOpt(),
    root_optimizer = StandardRootOpt(10),
    models_optimizer::ModelsUpdate = StandardModelsUpdate(),
) = StandardUpdate(
    nni,
    branchlength,
    root,
    models,
    refresh,
    argmax,
    branchlength_optimizer,
    root_optimizer,
    models_optimizer,
)

"""
    refresh!(tree::FelNode, models; partition_list = 1:length(tree.message))

Run `felsenstein!` and `felsenstein_down!` on `tree` with `models` to refresh messages.
"""
function refresh!(tree::FelNode, models; partition_list = 1:length(tree.message))
    felsenstein!(tree, models, partition_list = partition_list)
    felsenstein_down!(tree, models, partition_list = partition_list)
end

refresh!(update::StandardUpdate, tree::FelNode, models; partition_list = 1:length(tree.message)) = update.refresh && refresh!(tree, models, partition_list = partition_list)

function (update::StandardUpdate)(
    tree::FelNode,
    models;
    partition_list = 1:length(tree.message),
)
    update.nni > 0 && refresh!(update, tree, models, partition_list = partition_list)
    for _ = 1:update.nni
        nni_update!(update.nni_selection, tree, models, partition_list = partition_list)
    end
    update.branchlength > 0 && refresh!(update, tree, models, partition_list = partition_list)
    for _ = 1:update.branchlength
        branchlength_update!(
            update.branchlength_modifier,
            tree,
            models,
            partition_list = partition_list,
        )
    end
    update.root > 0 && refresh!(update, tree, models, partition_list = partition_list)
    for _ = 1:update.root
        tree = root_update!(
            update.root_update,
            tree,
            models,
            partition_list = partition_list,
        )
    end
    update.models > 0 && refresh!(update, tree, models, partition_list = partition_list)
    for _ = 1:update.models
        models = update.models_update(tree, models, partition_list = partition_list)
    end
    return tree, models
end

function collapse_models(u::StandardUpdate, models)
    return collapse_models(u.models_update, models)
end