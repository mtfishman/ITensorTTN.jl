using ITensorTTN
using Documenter

DocMeta.setdocmeta!(ITensorTTN, :DocTestSetup, :(using ITensorTTN); recursive=true)

makedocs(;
    modules=[ITensorTTN],
    authors="Matthew Fishman <mfishman@flatironinstitute.org> and contributors",
    repo="https://github.com/mtfishman/ITensorTTN.jl/blob/{commit}{path}#{line}",
    sitename="ITensorTTN.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mtfishman.github.io/ITensorTTN.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mtfishman/ITensorTTN.jl",
    devbranch="main",
)
