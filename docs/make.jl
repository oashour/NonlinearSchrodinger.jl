using NLSS
using Documenter

makedocs(;
    modules=[NLSS],
    authors="Omar Ashour",
    repo="https://github.com/oashour/NLSS.jl/blob/{commit}{path}#L{line}",
    sitename="NLSS.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://oashour.github.io/NLSS.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/oashour/NLSS.jl",
)
