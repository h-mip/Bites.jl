using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(path=joinpath(@__DIR__, ".."))
Pkg.instantiate()

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
        canonical="https://h-mip.github.io/Bites.jl/stable/",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

# Create a root index.html that redirects to stable
if get(ENV, "CI", "false") == "true"
    open(joinpath(@__DIR__, "build", "index.html"), "w") do io
        write(io, """
        <!DOCTYPE html>
        <html>
        <head>
            <meta http-equiv="refresh" content="0; url=stable/" />
            <meta name="robots" content="noindex">
            <link rel="canonical" href="stable/" />
        </head>
        <body>
            <p>Redirecting to <a href="stable/">stable documentation</a>...</p>
        </body>
        </html>
        """)
    end
end

if get(ENV, "CI", "false") == "true"
    deploydocs(;
        target = "build",
        repo="github.com/h-mip/Bites.jl",
        devbranch="main",
        devurl = "dev",
        versions = ["stable" => "v^", "v#.#"],
    )
end
