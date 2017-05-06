__precompile__()

module DiffEqNoiseProcess

using DataStructures, ResettableStacks, DiffEqBase, RecipesBase,
      RecursiveArrayTools

import DiffEqBase: isinplace, solve

include("types.jl")
include("interface.jl")
include("wiener.jl")
include("solve.jl")
include("geometric_bm.jl")
include("rswm.jl")
include("noise_wrapper_interface.jl")
include("recipes.jl")
include("correlated_noisefunc.jl")

export RSWM

export NoiseProcess, adaptive_alg

export WienerProcess, WienerProcess!, GeometricBrownianMotionProcess, GeometricBrownianMotionProcess!

export NoiseWrapper

export accept_step!, reject_step!, calculate_step!

export CorrelatedWienerProcess, CorrelatedWienerProcess!

end # module
