using DiffEqNoiseProcess, DiffEqBase, DiffEqMonteCarlo,
      Base.Test, DataStructures

W = BrownianBridge(0.0,1.0,0.0,1.0,0.0,0.0)
prob = NoiseProblem(W,(0.0,1.0))
monte_prob = MonteCarloProblem(prob)
@time sol = solve(monte_prob,dt=0.1,num_monte=10000)

# Spot check the mean and the variance
q =  0.4
@test ≈(timestep_mean(sol,5),q,atol=1e-2)
@test ≈(timestep_meanvar(sol,5)[2],(1-q)*q,atol=1e-2)
@test ≈(timestep_mean(sol,11)[1],1.0,atol=1e-16)
@test ≈(timestep_meanvar(sol,11)[2],0.0,atol=1e-16)

const μ = 1.2
const σ = 2.2
W = GeometricBrownianBridge(μ,σ,0.0,1.0,0.0,1.0,0.0,0.0)
prob = NoiseProblem(W,(0.0,1.0))
monte_prob = MonteCarloProblem(prob)
@time sol = solve(monte_prob,dt=0.1,num_monte=100)
