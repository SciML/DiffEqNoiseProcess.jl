using Documenter, DiffEqNoiseProcess


include("pages.jl")

makedocs(sitename = "DiffEqNoiseProcess.jl",
    authors = "Chris Rackauckas",
    modules = [DiffEqNoiseProcess],
    clean = true, doctest = false, linkcheck = true,
    warnonly = [:docs_block, :missing_docs],
    format = Documenter.HTML(assets = ["assets/favicon.ico"],
        canonical = "https://docs.sciml.ai/DiffEqNoiseProcess/stable/"),
    pages = pages)

deploydocs(repo = "github.com/SciML/DiffEqNoiseProcess.jl.git";
    push_preview = true)
