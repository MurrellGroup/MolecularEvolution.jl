using MolecularEvolution
using Documenter, Literate
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
    repo = "https://github.com/MurrellGroup/MolecularEvolution.jl/blob/{commit}{path}#{line}",
    sitename = "MolecularEvolution.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://MurrellGroup.github.io/MolecularEvolution.jl",
        edit_link = "main",
        assets = ["assets/favicon.ico"],
    ),
    pages = [
        "Home" => "index.md",
        "framework.md",
        "examples.md",
        "IO.md",
        "models.md",
        "simulation.md",
        "optimization.md",
        "ancestors.md",
        "generated/viz.md",
        "generated/update.md",
    ],
)

deploydocs(; repo = "github.com/MurrellGroup/MolecularEvolution.jl", devbranch = "main")
