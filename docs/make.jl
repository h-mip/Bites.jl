using Bites
using Documenter

DocMeta.setdocmeta!(Bites, :DocTestSetup, :(using Bites); recursive=true)

makedocs(;
    modules=[Bites],
    authors="John Palmer <johnrbpalmer@gmail.com> and contributors",
    repo="https://github.com/h-mip/Bites.jl/blob/{commit}{path}#{line}",
    sitename="Bites.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://h-mip.github.io/Bites.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/h-mip/Bites.jl",
    devbranch="main",
)
