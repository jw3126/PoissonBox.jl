using PoissonBox
using Documenter

DocMeta.setdocmeta!(PoissonBox, :DocTestSetup, :(using PoissonBox); recursive=true)

makedocs(;
    modules=[PoissonBox],
    authors="Jan Weidner <jw3126@gmail.com> and contributors",
    repo="https://github.com/jw3126/PoissonBox.jl/blob/{commit}{path}#{line}",
    sitename="PoissonBox.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jw3126.github.io/PoissonBox.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jw3126/PoissonBox.jl",
    devbranch="main",
)
