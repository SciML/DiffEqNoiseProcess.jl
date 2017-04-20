__precompile__()

module DiffEqNoiseProcess

using DataStructures, ResettableStacks, DiffEqBase

import DiffEqBase: isinplace

include("types.jl")
include("interface.jl")
include("wiener.jl")
include("rswm.jl")

export RSWM

export NoiseProcess, adaptive_alg, WienerProcess, WienerProcess!

export accept_step!, reject_step!, calculate_step!

end # module
