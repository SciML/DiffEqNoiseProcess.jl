@testset "OU" begin

using DiffEqNoiseProcess, DiffEqBase, DiffEqMonteCarlo, Test, Statistics
Θ = 1.0
μ = 1.2
σ = 0.3

W = OrnsteinUhlenbeckProcess(Θ,μ,σ,0.0,2.0,1.0)
prob = NoiseProblem(W,(0.0,1.0))
sol = solve(prob;dt=0.1)

t = 1.0
u0 = 2.0
expected_mean = μ + (u0 - μ)*exp(-Θ*t)
expected_variance = (1-exp(-2Θ*t))*σ^2/(2Θ)
monte_prob = MonteCarloProblem(prob;output_func = (sol,i)-> (sol[end],false))
sol = solve(monte_prob;dt=0.1,num_monte=10000)
@test abs(mean(sol) - expected_mean) < 0.04
@test abs(var(sol) - expected_variance) < 0.04

W = OrnsteinUhlenbeckProcess!(Θ,μ,σ,0.0,ones(2),ones(2))
prob = NoiseProblem(W,(0.0,1.0))
sol = solve(prob;dt=0.1)

end
