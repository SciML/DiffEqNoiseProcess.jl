using StochasticDiffEq, DiffEqNoiseProcess, Test, LinearAlgebra
using Random, Random123
seed = 100;
Random.seed!(seed);
@testset "Two noise processes for different m" begin
    for m in 1:4
        rand_prototype = zeros(m)
        rand_prototype2 = similar(rand_prototype, Int(m * (m - 1) / 2))
        rand_prototype2 .= false

        rng_base = Random123.Threefry4x()
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
