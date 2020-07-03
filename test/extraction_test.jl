@testset "Noise Extraction Test" begin
  using Random, DiffEqNoiseProcess, DiffEqBase, Test
  seed = 100
  Random.seed!(seed)

  tstart = 0.0
  tend = 1.0
  dt = 0.1
  trange = (tstart, tend)
  t = tstart:dt:tend
  tarray = collect(t)

  W = WienerProcess(0.0,0.0,0.0)
  prob = NoiseProblem(W,trange)
  sol = solve(prob, dt=dt)

  _sol = deepcopy(sol)
  sol.save_everystep = false
  for i in 1:1:length(tarray)
    t = tarray[i]
    sol(t)
  end
  @test length(_sol) == length(sol.W) == length(tarray)
  @test _sol == sol

  Random.seed!(seed)
  W2 = WienerProcess!(0.0,[0.0],[0.0])
  prob2 = NoiseProblem(W2,trange)
  sol2 = solve(prob2, dt=dt)

  sol2.save_everystep = false
  for i in 1:1:length(tarray)
    t = tarray[i]
    sol2(t)
  end

  @test length(sol2.W) == length(tarray)
  @test minimum(_sol .== sol2)
end
