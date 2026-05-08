using Test
using DiffEqNoiseProcess
using ChainRulesCore: ChainRulesCore, NoTangent, ProjectTo

@testset "ChainRulesCore: ProjectTo(::AbstractNoiseProcess) → NoTangent" begin
    ts = collect(0.0:0.1:1.0)

    # NoiseGrid with vector noise (the SDE2 / SciMLSensitivity#1432 trigger).
    Ws_vec = [randn(2) for _ in ts]
    grid_vec = NoiseGrid(ts, Ws_vec)
    @test ChainRulesCore.ProjectTo(grid_vec) isa ProjectTo{NoTangent}

    # NoiseGrid with scalar noise.
    Ws_scalar = randn(length(ts))
    grid_scalar = NoiseGrid(ts, Ws_scalar)
    @test ChainRulesCore.ProjectTo(grid_scalar) isa ProjectTo{NoTangent}

    # NoiseProcess (the canonical Wiener wrapper).
    np = WienerProcess(0.0, 0.0, nothing)
    @test ChainRulesCore.ProjectTo(np) isa ProjectTo{NoTangent}

    # Applying the override to a tangent yields NoTangent.
    proj = ChainRulesCore.ProjectTo(grid_vec)
    @test proj(rand(11)) isa NoTangent
    @test proj(NoTangent()) isa NoTangent
end
