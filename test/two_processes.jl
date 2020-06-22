using StochasticDiffEq, DiffEqNoiseProcess, Test, LinearAlgebra
using Random, RandomNumbers
seed=100; Random.seed!(seed)
@testset "Two noise processes for different m" begin
  for m=1:4
    rand_prototype = zeros(m)
    rand_prototype2 = similar(rand_prototype,Int(m*(m-1)/2))
    rand_prototype2 .= false
    Random.seed!(seed)
    W = WienerProcess!(0.0,rand_prototype,rand_prototype2,
                     save_everystep=true,
                     rng = Xorshifts.Xoroshiro128Plus(seed))
    prob = NoiseProblem(W,(0.0,1.0))
    sol = solve(prob;dt=0.2)

    Random.seed!(seed)
    Woop = WienerProcess(0.0,rand_prototype,rand_prototype2,
                     save_everystep=true,
                     rng = Xorshifts.Xoroshiro128Plus(seed))
    proboop = NoiseProblem(Woop,(0.0,1.0))
    soloop = solve(proboop;dt=0.2)

    @test norm(soloop.W - sol.W) == 0.0
    @test norm(soloop.Z - sol.Z) == 0.0

  end
end
