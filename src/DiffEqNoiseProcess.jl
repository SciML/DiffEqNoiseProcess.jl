module DiffEqNoiseProcess

using ResettableStacks, DiffEqBase, RecipesBase
using RecursiveArrayTools, StaticArraysCore, Random, Statistics
using LinearAlgebra

import RandomNumbers: Xorshifts

import RandomNumbers, Random123

import DiffEqBase: isinplace,
    solve, AbstractNoiseProcess,
    DEIntegrator, AbstractNoiseProblem

import PoissonRandom, Distributions

import QuadGK, Optim

import GPUArraysCore

using Markdown

using DiffEqBase: @..

include("types.jl")
include("copy_noise_types.jl")
include("wiener.jl")
include("solve.jl")
include("geometric_bm.jl")
include("compoundpoisson.jl")
include("ornstein_uhlenbeck.jl")
include("rswm.jl")
include("bridges.jl")
include("noise_interfaces/simple_noise_process_interface.jl")
include("noise_interfaces/noise_process_interface.jl")
include("noise_interfaces/noise_function_interface.jl")
include("noise_interfaces/virtual_brownian_tree_interface.jl")
include("noise_interfaces/noise_grid_interface.jl")
include("noise_interfaces/noise_approximation_interface.jl")
include("noise_interfaces/noise_wrapper_interface.jl")
include("noise_interfaces/box_wedge_tail_interface.jl")
include("noise_interfaces/noise_transport_interface.jl")
include("noise_interfaces/common.jl")
include("correlated_noisefunc.jl")
include("pCN.jl")

import Requires
@static if !isdefined(Base, :get_extension)
    function __init__()
        Requires.@require ReverseDiff="37e2e3b7-166d-5795-8a7a-e32c996b4267" begin
            include("../ext/DiffEqNoiseProcessReverseDiffExt.jl")
        end
    end
end

export RSWM

export NoiseProcess, SimpleNoiseProcess, adaptive_alg

export WienerProcess, WienerProcess!, SimpleWienerProcess, SimpleWienerProcess!

export RealWienerProcess, RealWienerProcess!

export BrownianBridge, BrownianBridge!

export GeometricBrownianMotionProcess, GeometricBrownianMotionProcess!

export CompoundPoissonProcess, CompoundPoissonProcess!

export GeometricBrownianBridge, GeometricBrownianBridge!

export OrnsteinUhlenbeckBridge, OrnsteinUhlenbeckBridge!

export CompoundPoissonBridge, CompoundPoissonBridge!

export OrnsteinUhlenbeckProcess, OrnsteinUhlenbeckProcess!

export NoiseWrapper, NoiseFunction, NoiseGrid, NoiseApproximation, NoiseTransport

export VirtualBrownianTree, VirtualBrownianTree!

export BoxWedgeTail, BoxWedgeTail!

export accept_step!, reject_step!, calculate_step!, setup_next_step!, save_noise!

export CorrelatedWienerProcess, CorrelatedWienerProcess!

export pCN, pCN!

end # module
