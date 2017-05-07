__precompile__()

module DiffEqNoiseProcess

using DataStructures, ResettableStacks, DiffEqBase, RecipesBase,
      RecursiveArrayTools

import DiffEqBase: isinplace, solve

include("types.jl")
include("wiener.jl")
include("solve.jl")
include("geometric_bm.jl")
include("ornstein_uhlenbeck.jl")
include("rswm.jl")
include("noise_interfaces/noise_process_interface.jl")
include("noise_interfaces/noise_function_interface.jl")
include("noise_interfaces/noise_grid_interface.jl")
include("noise_interfaces/noise_approximation_interface.jl")
include("noise_interfaces/noise_wrapper_interface.jl")
include("recipes.jl")
include("correlated_noisefunc.jl")

export RSWM

export NoiseProcess, adaptive_alg

export WienerProcess, WienerProcess!

export GeometricBrownianMotionProcess, GeometricBrownianMotionProcess!

export OrnsteinUhlenbeckProcess, OrnsteinUhlenbeckProcess!

export NoiseWrapper, NoiseFunction, NoiseGrid, NoiseApproximation

export accept_step!, reject_step!, calculate_step!

export CorrelatedWienerProcess, CorrelatedWienerProcess!

end # module
