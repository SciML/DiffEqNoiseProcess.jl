module DiffEqNoiseProcess

using ResettableStacks: ResettableStacks
using DiffEqBase: DiffEqBase
using RecipesBase: RecipesBase
using RecursiveArrayTools: RecursiveArrayTools, recursivecopy
using StaticArraysCore: StaticArraysCore, SArray
using Random: Random, AbstractRNG, randn!
using Statistics: Statistics
using LinearAlgebra: LinearAlgebra, Diagonal, mul!, svd

import Random123

import DiffEqBase: isinplace, AbstractNoiseProcess,
    DEIntegrator, AbstractNoiseProblem

import SciMLBase
import SciMLBase: add_tstop!, reinit!

import CommonSolve: step!

import PoissonRandom, Distributions

import QuadGK

import GPUArraysCore

using Markdown: Markdown, @doc_str

using DiffEqBase: @..

using Base: deleteat!, convert, copyto!

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
