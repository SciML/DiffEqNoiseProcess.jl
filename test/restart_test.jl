@testset "Restart" begin

  using Test, LinearAlgebra
  using StochasticDiffEq, DiffEqNoiseProcess
  using Random

  u₀ = [0.5]
  tstart = 0.0
  tend = 1.0
  dt = 0.005
  trange = (tstart, tend)

  f!(du,u,p,t) = du.=1.01*u
  σ!(du,u,p,t) = du.=0.87*u

  dt1 = tend/1e3
  seed = 100
  Random.seed!(seed)
  prob = SDEProblem(f!,σ!,u₀,trange)
  sol = solve(prob,EulerHeun(),dt=dt1,adaptive=false, save_noise=true, saveat=collect(trange[1]:dt1:trange[2]))

  # choose a random interval
  interval = (0.35, 0.87)

  idx1 = searchsortedfirst(sol.t, interval[1])
  idx2 = searchsortedfirst(sol.t, interval[2])

  forwardnoise = DiffEqNoiseProcess.NoiseGrid(sol.t[idx1:idx2], sol.W.W[idx1:idx2])
  forwardnoise2 = DiffEqNoiseProcess.NoiseWrapper(sol.W, indx=idx1)

  cpsol = solve(remake(prob, tspan=interval, u0=sol(interval[1]), noise=forwardnoise), sol.alg, save_noise=false; dt=dt1)
  cpsol2 = solve(remake(prob, tspan=interval, u0=sol(interval[1]), noise=forwardnoise2), sol.alg, save_noise=false; dt=dt1)

  sola = vcat(sol.u[idx1:idx2] ...)
  cpsola = vcat(cpsol.u ...)
  cpsol2a = vcat(cpsol2.u ...)

  @test isapprox(sola, cpsola, atol=1e-10 )
  @test isapprox(sola, cpsol2a, atol=1e-10 )

end
