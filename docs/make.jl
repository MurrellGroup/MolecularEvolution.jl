using MolecularEvolution
using Documenter

DocMeta.setdocmeta!(MolecularEvolution, :DocTestSetup, :(using MolecularEvolution); recursive=true)

makedocs(;
    modules=[MolecularEvolution],
    authors="Ben Murrell <murrellb@gmail.com> and contributors",
    repo="https://github.com/murrellb/MolecularEvolution.jl/blob/{commit}{path}#{line}",
    sitename="MolecularEvolution.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://murrellb.github.io/MolecularEvolution.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/murrellb/MolecularEvolution.jl",
    devbranch="main",
)
