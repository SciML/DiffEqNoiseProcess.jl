using Documenter, DiffEqNoiseProcess

makedocs(
    sitename = "DiffEqNoiseProcess.jl Test",
    modules = [DiffEqNoiseProcess],
    format = Documenter.HTML(
        prettyurls = false
    ),
    pages = [
        "API" => "api/interface.md"
    ],
    warnonly = true
)