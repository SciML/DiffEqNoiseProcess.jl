using Documenter, DiffEqNoiseProcess

include("pages.jl")

makedocs(sitename = "DiffEqNoiseProcess.jl",
         authors = "Chris Rackauckas",
         modules = [DiffEqNoiseProcess],
         clean = true, doctest = false,
         format = Documenter.HTML(analytics = "UA-90474609-3",
                                  assets = ["assets/favicon.ico"],
                                  canonical = "https://noise.sciml.ai/stable/"),
         pages = pages)

deploydocs(repo = "github.com/SciML/DiffEqNoiseProcess.jl.git";
           push_preview = true)
