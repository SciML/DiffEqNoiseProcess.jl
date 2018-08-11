@testset "Brownian Bridge" begin

using DiffEqNoiseProcess, DiffEqBase, DiffEqMonteCarlo,
      Test, DataStructures, Random

Random.seed!(100)
W = BrownianBridge(0.0,1.0,0.0,1.0,0.0,0.0)
prob = NoiseProblem(W,(0.0,1.0))
monte_prob = MonteCarloProblem(prob)
@time sol = solve(monte_prob,dt=0.1,num_monte=10000)

# Spot check the mean and the variance
qs = 0:0.1:1
for i in 2:10
  q = qs[i]
  @test ≈(timestep_mean(sol,i),q,atol=1e-2)
  @test ≈(timestep_meanvar(sol,i)[2],(1-q)*q,atol=1e-2)
end
@test ≈(timestep_mean(sol,1)[1],0.0,atol=1e-16)
@test ≈(timestep_meanvar(sol,1)[2],0.0,atol=1e-16)
@test ≈(timestep_mean(sol,11)[1],1.0,atol=1e-16)
@test ≈(timestep_meanvar(sol,11)[2],0.0,atol=1e-16)


μ = 1.2
σ = 2.2
W = GeometricBrownianBridge(μ,σ,0.0,1.0,0.0,1.0,0.0,0.0)
prob = NoiseProblem(W,(0.0,1.0))
monte_prob = MonteCarloProblem(prob)
@time sol = solve(monte_prob,dt=0.1,num_monte=100)

end
