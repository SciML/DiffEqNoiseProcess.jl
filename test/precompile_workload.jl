using DiffEqNoiseProcess
using Random
using SciMLBase: AbstractNoiseProcess
using Test

@testset "Precompile workload" begin
    process = WienerProcess(0.0, 0.0; rng = Xoshiro(0))
    @test process isa AbstractNoiseProcess

    calculate_step!(process, 0.1, nothing, nothing)
    accept_step!(process, 0.1, nothing, nothing, false)
    value, auxiliary = process(0.05)
    @test value isa Float64
    @test auxiliary === nothing

    reject_step!(process, 0.05, nothing, nothing)
    @test accept_step!(process, 0.05, nothing, nothing, false) === nothing

    in_place = WienerProcess!(0.0, zeros(2); rng = Xoshiro(0))
    calculate_step!(in_place, 0.1, nothing, nothing)
    accept_step!(in_place, 0.1, nothing, nothing, false)
    value = in_place(0.05)[1]
    @test value isa Vector{Float64}
    @test length(value) == 2
end
