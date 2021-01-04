using PRS
using Documenter

makedocs(;
    modules=[PRS],
    authors="Florian Zink <zink.florian@gmail.com> and contributors",
    repo="https://github.com/mfz/PRS.jl/blob/{commit}{path}#L{line}",
    sitename="PRS.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mfz.github.io/PRS.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mfz/PRS.jl",
)
