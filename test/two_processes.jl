using DiffEqBase, DiffEqNoiseProcess, Test, LinearAlgebra
using Random
seed = 100;
Random.seed!(seed);
@testset "Two noise processes for different m" begin
    for m in 1:4
        rand_prototype = zeros(m)
        rand_prototype2 = similar(rand_prototype, Int(m * (m - 1) / 2))
        rand_prototype2 .= false

        rng_base = Xoshiro()
        Random.seed!(rng_base)
        rng1 = copy(rng_base)
        rng2 = copy(rng_base)

        W = WienerProcess!(
            0.0, copy(rand_prototype), copy(rand_prototype2),
            save_everystep = true,
            rng = rng1, reseed = false
        )
        prob = NoiseProblem(W, (0.0, 1.0))
        sol = solve(prob; dt = 0.2)

        Woop = WienerProcess(
            0.0, copy(rand_prototype), copy(rand_prototype2),
            save_everystep = true,
            rng = rng2, reseed = false
        )
        proboop = NoiseProblem(Woop, (0.0, 1.0))
        soloop = solve(proboop; dt = 0.2)

        @test norm(soloop.W - sol.W) == 0.0
        @test norm(soloop.Z - sol.Z) == 0.0
    end
end

@testset "Default RNG produces independent noise processes" begin
    # reseed = false so that solve doesn't re-randomize the RNG,
    # ensuring we test that construction-time RNG states differ.
    W1 = WienerProcess(0.0, 0.0, 0.0; reseed = false)
    W2 = WienerProcess(0.0, 0.0, 0.0; reseed = false)
    sol1 = solve(NoiseProblem(W1, (0.0, 1.0)); dt = 0.1)
    sol2 = solve(NoiseProblem(W2, (0.0, 1.0)); dt = 0.1)
    @test sol1.W != sol2.W

    W1! = WienerProcess!(0.0, zeros(3), zeros(3); reseed = false)
    W2! = WienerProcess!(0.0, zeros(3), zeros(3); reseed = false)
    sol1! = solve(NoiseProblem(W1!, (0.0, 1.0)); dt = 0.1)
    sol2! = solve(NoiseProblem(W2!, (0.0, 1.0)); dt = 0.1)
    @test sol1!.W != sol2!.W
end
