struct StandardModelsUpdate <: ModelsUpdate end

(::StandardModelsUpdate)(tree::FelNode, models; partition_list = 1:length(tree.message)) =
    models

collapse_models(::ModelsUpdate, models) = models
function collapse_models(::StandardModelsUpdate, models) end