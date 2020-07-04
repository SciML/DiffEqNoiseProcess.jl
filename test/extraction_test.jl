@testset "Noise Extraction Test" begin
  using Random, DiffEqNoiseProcess, DiffEqBase, Test
  seed = 100
  Random.seed!(seed)

  tstart = 0.0
  tend = 1.0
  dt = 0.001
  trange = (tstart, tend)

  Random.seed!(seed)
  W = WienerProcess(0.0,0.0,0.0)
  prob = NoiseProblem(W,trange)
  sol = solve(prob, dt=dt)

  _sol = deepcopy(sol)
  _sol.save_everystep = false
  sol.save_everystep = false
  tarray = [tstart]
  for i = 1:1000
    t = tarray[i]
    sol(t)
    push!(tarray,t+dt)
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


  W3 = NoiseGrid(reverse(_sol.t),reverse(_sol.W))
  prob3 = NoiseProblem(W3,reverse(trange))
  sol3 = solve(prob3, dt=-dt)

  @test minimum(reverse(sol3.W) .== _sol)

  W4 = NoiseWrapper(_sol, reverse=true)
  for i in 1:1:length(tarray)
    t = tarray[end-i+1]
    tmp = DiffEqNoiseProcess.interpolate!(W4,nothing,nothing,t)
    @test _sol.W[end-i+1] ≈ tmp[1]
    @test _sol.Z[end-i+1] ≈ tmp[2]
  end

end
