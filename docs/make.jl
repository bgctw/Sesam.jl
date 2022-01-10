using Sesam
using Documenter

DocMeta.setdocmeta!(Sesam, :DocTestSetup, :(using Sesam); recursive=true)

makedocs(;
    modules=[Sesam],
    authors="Thomas Wutzler <twutz@bgc-jena.mpg.de> and contributors",
    repo="https://github.com/bgctw/Sesam.jl/blob/{commit}{path}#{line}",
    sitename="Sesam.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://bgctw.github.io/Sesam.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/bgctw/Sesam.jl",
    devbranch="main",
)
