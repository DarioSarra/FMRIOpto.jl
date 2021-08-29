using Documenter, FMRIOpto

makedocs(;
    modules=[FMRIOpto],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/DarioSarra/FMRIOpto.jl/blob/{commit}{path}#L{line}",
    sitename="FMRIOpto.jl",
    authors="DarioSarra",
    assets=String[],
)

deploydocs(;
    repo="github.com/DarioSarra/FMRIOpto.jl",
)
