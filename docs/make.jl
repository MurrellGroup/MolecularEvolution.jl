using MolecularEvolution
using Documenter, DocumenterVitepress, Literate
using Phylo, Distributions
using Plots
using Compose, Cairo, Fontconfig
using FASTX

LITERATE_INPUT = joinpath(@__DIR__, "..", "examples")
LITERATE_OUTPUT = OUTPUT = joinpath(@__DIR__, "src/generated")

for (root, _, files) ∈ walkdir(LITERATE_INPUT), file ∈ files
    @show root, files
    # ignore non julia files
    splitext(file)[2] == ".jl" || continue
    # full path to a literate script
    ipath = joinpath(root, file)
    # generated output path
    opath = splitdir(replace(ipath, LITERATE_INPUT=>LITERATE_OUTPUT))[1]
    # generate the markdown file calling Literate
    Literate.markdown(ipath, opath)
end

DocMeta.setdocmeta!(
    MolecularEvolution,
    :DocTestSetup,
    :(using MolecularEvolution);
    recursive = true,
)

makedocs(;
    modules = [MolecularEvolution],
    authors = "Ben Murrell <benjamin.murrell@ki.se> and contributors",
    sitename = "MolecularEvolution.jl",
    format = DocumenterVitepress.MarkdownVitepress(;
        repo = "https://github.com/MurrellGroup/MolecularEvolution.jl",
        #=prettyurls = get(ENV, "CI", "false") == "true",=#
        deploy_url = "https://MurrellGroup.github.io/MolecularEvolution.jl",
        devbranch = "main",
        devurl = "dev",
        assets = ["assets/favicon.ico"],
    #=for local hosting
    
        md_output_path = ".",
        build_vitepress = false
    ),
    clean = false,
    =#
    ),
    pages = [
        "Core Concepts" => [
            "intro.md",
            "framework.md",
            "models.md",
        ],
        "Methods & Algorithms" => [
            "simulation.md",
            "optimization.md",
            "ancestors.md",
        ],
        "Extensions & Utilities" => [
            "IO.md",
            "generated/viz.md",
            "generated/update.md",
        ],
        "examples.md",
        "api.md"
    ],
)

deploydocs(;
    repo = "github.com/MurrellGroup/MolecularEvolution.jl",
    target = "build", # this is where Vitepress stores its output
    devbranch = "main",
    branch = "gh-pages",
    push_preview = true,
)
