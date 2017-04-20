module DiffEqNoiseProcess

using DiffEqBase, DataStructures, ResettableStacks, StochasticDiffEq

import DiffEqBase: isinplace

include("types.jl")
include("interface.jl")
include("wiener.jl")

export NoiseProcess, adaptive_alg, WienerProcess, WienerProcess!

export accept_step!, reject_step!, calculate_step!

end # module
