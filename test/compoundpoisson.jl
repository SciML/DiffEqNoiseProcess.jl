@testset "CompoundPoissonProcess" begin

using DiffEqNoiseProcess, DiffEqBase, Test, Statistics

r = 2
rate(u,p,t) = r
W = CompoundPoissonProcess(rate,0.0,1.0)
dt = 0.1
calculate_step!(W,dt,nothing,nothing)
for i in 1:10
  accept_step!(W,dt,nothing,nothing)
end

W = CompoundPoissonProcess(rate,0.0,1.0)
prob = NoiseProblem(W,(0.0,1.0))
sol = solve(prob;dt=0.1)

W = CompoundPoissonProcess(rate,0//1,1.0)
prob_rational = NoiseProblem(W,(0,1))
dt = 1//10
sol = solve(prob_rational;dt=dt)

dt = dt/100
for t in dt:dt:1-dt
  sol(t)
end

t = 1.0
u0 = 1.0
expected_mean = 1 + t*r
expected_variance = t*r
ensemble_prob = EnsembleProblem(prob;output_func = (sol,i)-> (sol[end],false))
sol = solve(ensemble_prob;dt=0.1,trajectories=40000)
@test abs(mean(sol) - expected_mean) < 0.4
@test abs(var(sol) - expected_variance) < 0.4 # Converges slowly

rate(du,u,p,t) = (@. du = r)
W = CompoundPoissonProcess!(rate,0.0,ones(2))
prob = NoiseProblem(W,(0.0,1.0))
sol = solve(prob;dt=0.1)

ensemble_prob = EnsembleProblem(prob;output_func = (sol,i)-> (sol[end],false))
sol = solve(ensemble_prob;dt=0.1,trajectories=40000)
@test abs(mean(sol) - expected_mean) < 0.4
@test abs(var(sol) - expected_variance) < 0.4 # Converges slowly
end
