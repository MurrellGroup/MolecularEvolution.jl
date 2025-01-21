struct AbstractUpdater
    nni::Bool
    branchlength::Bool
    root::Bool
    nni_selection::Function
    branchlength_sampler::UnivariateModifier
    root_sampler::Any
    AbstractUpdater(
        nni::Bool,
        branchlength::Bool,
        root::Bool,
        nni_selection::Function,
        branchlength_sampler::UnivariateModifier,
        root_sampler::Any,
    ) = new(nni, branchlength, root, nni_selection, branchlength_sampler, root_sampler)
end

Bayes(;nni=true, branchlength=true, root=true, nni_selection=softmax_selection, branchlength_sampler=BranchlengthSampler(Normal(0,2), Normal(-1,1)), root_sampler=RootPositionSampler()) = AbstractUpdater(nni, branchlength, root, nni_selection, branchlength_sampler, root_sampler)

function update!(u::Updater, tree::FelNode, models)
    nni_update!(u, tree, models)
    branchlength_update!(u, tree, models)
    root_update!(u, tree, models)
    return tree
end
