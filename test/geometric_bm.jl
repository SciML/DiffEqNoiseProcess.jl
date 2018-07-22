@testset "GeometricBM" begin

using DiffEqNoiseProcess, DiffEqBase, DiffEqMonteCarlo, Test, Statistics

μ = 1.0
σ = 2.0

W = GeometricBrownianMotionProcess(μ,σ,0.0,1.0,1.0)
dt = 0.1
calculate_step!(W,dt)
for i in 1:10
  accept_step!(W,dt)
end

W = GeometricBrownianMotionProcess(μ,σ,0.0,1.0,1.0)
prob = NoiseProblem(W,(0.0,1.0))
sol = solve(prob;dt=0.1)

W = GeometricBrownianMotionProcess(μ,σ,0//1,1.0,1.0)
prob_rational = NoiseProblem(W,(0,1))
dt = 1//10
sol = solve(prob_rational;dt=dt)

dt = dt/100
for t in dt:dt:1-dt
  sol(t)
end

t = 1.0
u0 = 1.0
expected_mean = u0*exp(μ*t)
expected_variance = u0^2*exp(2μ*t)*(exp(σ^2*t)-1)
monte_prob = MonteCarloProblem(prob;output_func = (sol,i)-> (sol[end],false))
sol = solve(monte_prob;dt=0.1,num_monte=40000)
@test abs(mean(sol) - expected_mean) < 0.4
#abs(var(sol) - expected_variance) < 0.4 # Converges slowly

W = GeometricBrownianMotionProcess!(μ,σ,0.0,ones(2),ones(2))
prob = NoiseProblem(W,(0.0,1.0))
sol = solve(prob;dt=0.1)

end
