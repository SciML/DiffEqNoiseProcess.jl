@testset "NoiseGrid" begin

using DiffEqNoiseProcess, DiffEqBase, Test

t = 0:0.001:1
grid = exp.(t)
W = NoiseGrid(t,grid)

dt = 0.1
calculate_step!(W,dt)

for i in 1:10
  accept_step!(W,dt)
end

W = NoiseGrid(t,grid)
prob = NoiseProblem(W,(0.0,1.0))
sol = solve(prob;dt=0.1)

@test !sol.step_setup
@test_throws ErrorException accept_step!(sol,dt)

t = 0:0.001:1
grid = [[exp.(ti) for i in 1:8] for ti in t]
W = NoiseGrid(t,grid)
prob = NoiseProblem(W,(0.0,1.0))
sol = solve(prob;dt=0.1)

dt = 0.001
t = 0:dt:1
brownian_values = cumsum([0;[sqrt(dt)*randn() for i in 1:length(t)-1]])
W = NoiseGrid(t,grid)

end
