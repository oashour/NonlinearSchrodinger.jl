using NLSS
using Documenter

makedocs(;
    modules=[NLSS],
    authors="Omar Ashour",
    repo="https://github.com/oashour/NLSS.jl/blob/{commit}{path}#L{line}",
    sitename="NLSS.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="http://oashour.github.io",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/oashour/NLSS.jl",
    branch = "gh-pages",
    devbranch = "master",
    devurl = "dev"
    #versions = ["stable" => "v^", "v#.#", devurl => devurl],
)
