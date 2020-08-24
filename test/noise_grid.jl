@testset "NoiseGrid" begin

  using DiffEqNoiseProcess, DiffEqBase, Test

  t = 0:0.001:1
  grid = exp.(t)
  W = NoiseGrid(t,grid)

  dt = 0.1
  calculate_step!(W,dt,nothing,nothing)

  for i in 1:10
    accept_step!(W,dt,nothing,nothing)
  end

  W = NoiseGrid(t,grid)
  prob = NoiseProblem(W,(0.0,1.0))
  sol = solve(prob;dt=0.1)

  @test !sol.step_setup
  @test_throws ErrorException accept_step!(sol,dt,nothing,nothing)

  t = 0:0.001:1
  grid = [[exp.(ti) for i in 1:8] for ti in t]
  W = NoiseGrid(t,grid)
  prob = NoiseProblem(W,(0.0,1.0))
  sol = solve(prob;dt=0.1)

  dt = 0.001
  t = 0:dt:1
  brownian_values = cumsum([0;[sqrt(dt)*randn() for i in 1:length(t)-1]])
  W = NoiseGrid(t,brownian_values)

  dt = 0.001
  t = 0:dt:1
  brownian_values2 = cumsum([[zeros(8)];[sqrt(dt)*randn(8) for i in 1:length(t)-1]])
  W = NoiseGrid(t,brownian_values2)
  prob = NoiseProblem(W,(0.0,1.0))
  sol = solve(prob;dt=0.1)

  dt = 1//1000
  t = 0:dt:1
  W = NoiseGrid(t,brownian_values)
  prob_rational = NoiseProblem(W,(0,1))
  sol = solve(prob_rational; dt=1//10)
end
