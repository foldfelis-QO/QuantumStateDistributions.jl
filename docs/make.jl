using QuantumStateDistributions
using Documenter

DocMeta.setdocmeta!(QuantumStateDistributions, :DocTestSetup, :(using QuantumStateDistributions); recursive=true)

makedocs(;
    modules=[QuantumStateDistributions],
    authors="JingYu Ning <foldfelis@gmail.com> and contributors",
    repo="https://github.com/foldfelis-QO/QuantumStateDistributions.jl/blob/{commit}{path}#{line}",
    sitename="QuantumStateDistributions.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://foldfelis-QO.github.io/QuantumStateDistributions.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/foldfelis-QO/QuantumStateDistributions.jl",
    devbranch="main",
)
