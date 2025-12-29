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
        canonical="https://h-mip.com/Bites.jl/stable/",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

# Create a root index.html landing page in the build directory
if get(ENV, "CI", "false") == "true"
    build_dir = joinpath(@__DIR__, "build")
    if isdir(build_dir)
        open(joinpath(build_dir, "index.html"), "w") do io
            write(io, """
            <!DOCTYPE html>
            <html lang="en">
            <head>
                <meta charset="UTF-8">
                <meta name="viewport" content="width=device-width, initial-scale=1.0">
                <title>Bites.jl Documentation</title>
                <style>
                    body {
                        font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, sans-serif;
                        max-width: 800px;
                        margin: 80px auto;
                        padding: 20px;
                        line-height: 1.6;
                    }
                    h1 {
                        color: #333;
                        border-bottom: 2px solid #4CAF50;
                        padding-bottom: 10px;
                    }
                    .version-list {
                        list-style: none;
                        padding: 0;
                    }
                    .version-list li {
                        margin: 15px 0;
                    }
                    .version-list a {
                        display: inline-block;
                        padding: 12px 24px;
                        background-color: #4CAF50;
                        color: white;
                        text-decoration: none;
                        border-radius: 4px;
                        transition: background-color 0.3s;
                        min-width: 200px;
                        text-align: center;
                    }
                    .version-list a:hover {
                        background-color: #45a049;
                    }
                    .description {
                        color: #666;
                        font-size: 0.9em;
                        margin-left: 10px;
                    }
                    .dev-version {
                        background-color: #2196F3;
                    }
                    .dev-version:hover {
                        background-color: #0b7dda;
                    }
                </style>
            </head>
            <body>
                <h1>Bites.jl Documentation</h1>
                <p>Welcome to the Bites.jl documentation. Please select a version:</p>
                
                <ul class="version-list">
                    <li>
                        <a href="stable/">📘 Stable Documentation</a>
                        <span class="description">Recommended for most users</span>
                    </li>
                    <li>
                        <a href="dev/" class="dev-version">🚧 Development Documentation</a>
                        <span class="description">Latest features from the main branch</span>
                    </li>
                </ul>
                
                <p style="margin-top: 40px; color: #666; font-size: 0.9em;">
                    For older versions, please visit the <a href="stable/" style="color: #4CAF50;">stable documentation</a> 
                    and use the version selector in the sidebar.
                </p>
            </body>
            </html>
            """)
        end
    else
        @warn "Build directory not found, skipping root index.html creation"
    end
end

# Deploy documentation
if get(ENV, "CI", "false") == "true"
    deploydocs(;
        repo="github.com/h-mip/Bites.jl",
        devbranch="main",
        devurl = "dev",
        versions = ["stable" => "v^", "v#.#"],
    )
end
