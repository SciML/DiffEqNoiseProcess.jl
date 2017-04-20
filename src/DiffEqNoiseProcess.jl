__precompile__()

module DiffEqNoiseProcess

using DataStructures, ResettableStacks, DiffEqBase, RecipesBase

import DiffEqBase: isinplace

include("types.jl")
include("interface.jl")
include("wiener.jl")
include("rswm.jl")
include("noise_wrapper_interface.jl")
include("recipes.jl")

export RSWM

export NoiseProcess, adaptive_alg, WienerProcess, WienerProcess!

export NoiseWrapper

export accept_step!, reject_step!, calculate_step!

end # module
