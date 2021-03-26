@testset "preconditioned Crank Nicolson tests" begin

using DiffEqNoiseProcess, Test, Random
using Statistics
using DiffEqBase
using DiffEqBase.EnsembleAnalysis

##
# Tests with
##

W = WienerProcess(0.0,0.0,0.0)

dt = 0.1
calculate_step!(W,dt,nothing,nothing)

for i in 1:10
  accept_step!(W,dt,nothing,nothing)
end

_W = deepcopy(W)

# test with ρ=1
W2 = pCN(_W, 1.0)
WWrapper = NoiseWrapper(_W)


# test source
@test W2.source.W == WWrapper.source.W
@test W2.source.Z == WWrapper.source.Z
@test W2.source.W == W.W
@test W2.source.Z == W.Z

# multi dimensional
W = WienerProcess(0.0,zeros(4),zeros(4))

dt = 0.1
calculate_step!(W,dt,nothing,nothing)

for i in 1:10
  accept_step!(W,dt,nothing,nothing)
end

_W = deepcopy(W)

W2 = pCN(_W, 1.0)
WWrapper = NoiseWrapper(_W)

# test source
@test W2.source.W == WWrapper.source.W
@test W2.source.Z == WWrapper.source.Z
@test W2.source.W == _W.W
@test W2.source.Z == _W.Z
@test W2.source.W == W.W
@test W2.source.Z == W.Z


# test ρ!=0 and ρ!=1
W3 = pCN(_W, 0.2)

# test source
@test W3.source.W != W.W
@test W3.source.Z == W.Z # no action on auxilary process

# inplace
W = WienerProcess!(0.0,zeros(4),zeros(4))

dt = 0.1
calculate_step!(W,dt,nothing,nothing)

for i in 1:10
  accept_step!(W,dt,nothing,nothing)
end

_W = deepcopy(W)

# test with ρ=0
W2 = pCN(_W, 0.0)
WWrapper = NoiseWrapper(_W)

# test source
@test W2.source.W != W.W
@test W2.source.Z == W.Z

# statistics test
W = WienerProcess(0.0,0.0,0.0)
prob = NoiseProblem(W,(0.0,1.0))
sol = solve(prob,dt=0.1)

function prob_func(prob,i,repeat)
  _sol = deepcopy(sol)
  Wtmp = pCN(_sol,1.0)
  remake(prob, noise=Wtmp)
end

ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)
@time sim = solve(ensemble_prob,dt=0.1,trajectories=100)

# Spot check the mean and the variance
qs = 0:0.1:1
for i in 1:11
  q = qs[i]
  @test ≈(timestep_mean(sim,i),sol.W[i],atol=1e-2)
  @test ≈(timestep_meanvar(sim,i)[2],zero(sol.W[i]),atol=1e-2)
end


# using Plots
# m_series = [timestep_mean(sim,i)  for i in 1:11]
# plot(m_series)
# plot!(sol.W)

end
