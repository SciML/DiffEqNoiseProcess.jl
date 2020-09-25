using StochasticDiffEq, DiffEqNoiseProcess, Test, Random
@testset "SDE Stratonovich Reversal Tests" begin
  Random.seed!(100)
  α=1.01
  β=0.87

  dt = 1e-3
  tspan = (0.0,1.0)
  u₀=1/2

  tarray =  collect(tspan[1]:dt:tspan[2])

  f!(du,u,p,t) = du .= α*u
  g!(du,u,p,t) = du .= β*u


  prob = SDEProblem(f!,g!,[u₀],tspan)
  sol = solve(prob,EulerHeun(),dt=dt,save_noise=true, adaptive=false)

  _sol = deepcopy(sol) # to make sure the plot is correct
  W1 = NoiseGrid(reverse!(_sol.t),reverse!(_sol.W.W))
  prob1 = SDEProblem(f!,g!,sol[end],reverse(tspan),noise=W1)
  sol1 = solve(prob1,EulerHeun(),dt=dt)

  _sol = deepcopy(sol)
  W2 = NoiseWrapper(_sol.W, reverse=true)
  prob2 = SDEProblem(f!,g!,sol[end],reverse(tspan),noise=W2)
  sol2 = solve(prob2,EulerHeun(),dt=dt, save_noise=false)

  # same time steps

  @test sol.u ≈ reverse(sol1.u) atol=1e-2
  @test sol.u ≈ reverse(sol2.u) atol=1e-2
  @test sol1.u ≈ sol2.u atol=1e-5

  # 1/10 size time steps

  _sol = deepcopy(sol) # to make sure the plot is correct
  W1 = NoiseGrid(reverse!(_sol.t),reverse!(_sol.W.W))
  prob1 = SDEProblem(f!,g!,sol[end],reverse(tspan),noise=W1)
  sol1 = solve(prob1,EulerHeun(),dt=0.1*dt)

  _sol = deepcopy(sol)
  W2 = NoiseWrapper(_sol.W, reverse=true)
  prob2 = SDEProblem(f!,g!,sol[end],reverse(tspan),noise=W2)
  sol2 = solve(prob2,EulerHeun(),dt=0.1*dt, save_noise=false)

  @test sol.u ≈ sol1(tarray).u atol=1e-2
  @test sol.u ≈ sol2(tarray).u atol=1e-2


  # diagonal noise

  prob = SDEProblem(f!,g!,[u₀,u₀],tspan)
  sol = solve(prob,EulerHeun(),dt=dt,save_noise=true)

  _sol = deepcopy(sol) # to make sure the plot is correct
  W1 = NoiseGrid(reverse!(_sol.t),reverse!(_sol.W.W))
  prob1 = SDEProblem(f!,g!,sol[end],reverse(tspan),noise=W1)
  sol1 = solve(prob1,EulerHeun(),dt=dt)

  _sol = deepcopy(sol)
  W2 = NoiseWrapper(_sol.W, reverse=true)
  prob2 = SDEProblem(f!,g!,sol[end],reverse(tspan),noise=W2)
  sol2 = solve(prob2,EulerHeun(),dt=dt)

  @test sol.u ≈ reverse(sol1.u) atol=5e-2
  @test sol.u ≈ reverse(sol2.u) atol=5e-2
  @test sol1.u ≈ sol2.u atol=1e-5


  # non-diagonal noise

  function gnd!(du,u,p,t)
    du[1,1] = 0.3u[1]
    du[1,2] = 0.6u[1]
    du[1,3] = 0.9u[1]
    du[1,4] = 0.12u[2]
    du[2,1] = 1.2u[1]
    du[2,2] = 0.2u[2]
    du[2,3] = 0.3u[2]
    du[2,4] = 1.8u[2]
  end
  prob = SDEProblem(f!,gnd!,[u₀,u₀],tspan,noise_rate_prototype=zeros(2,4))
  sol = solve(prob,EulerHeun(),dt=dt,save_noise=true)

  _sol = deepcopy(sol) # to make sure the plot is correct
  W1 = NoiseGrid(reverse!(_sol.t),reverse!(_sol.W.W))
  prob1 = SDEProblem(f!,gnd!,sol[end],reverse(tspan),noise=W1,noise_rate_prototype=zeros(2,4))
  sol1 = solve(prob1,EulerHeun(),dt=dt)

  _sol = deepcopy(sol)
  W2 = NoiseWrapper(_sol.W, reverse=true)
  prob2 = SDEProblem(f!,gnd!,sol[end],reverse(tspan),noise=W2, noise_rate_prototype=zeros(2,4))
  sol2 = solve(prob2,EulerHeun(),dt=dt)

  @test sol.u ≈ reverse(sol1.u) atol=5e-1
  @test sol.u ≈ reverse(sol2.u) atol=5e-1
  @test sol1.u ≈ sol2.u atol=1e-5

  ###
  ### OOP
  ###

  f(u,p,t) = α*u
  g(u,p,t) = β*u

  prob = SDEProblem(f,g,u₀,tspan)
  sol =solve(prob,EulerHeun(),dt=dt,save_noise=true)

  _sol = deepcopy(sol) # to make sure the plot is correct
  W1 = NoiseGrid(reverse!(_sol.t),reverse!(_sol.W.W))
  prob1 = SDEProblem(f,g,sol[end],reverse(tspan),noise=W1)
  sol1 = solve(prob1,EulerHeun(),dt=dt)

  _sol = deepcopy(sol)
  W2 = NoiseWrapper(_sol.W, reverse=true)
  prob2 = SDEProblem(f,g,sol[end],reverse(tspan),noise=W2)
  sol2 = solve(prob2,EulerHeun(),dt=dt)

  @test sol.u ≈ reverse(sol1.u) atol=2e-2
  @test sol.u ≈ reverse(sol2.u) atol=2e-2
  @test sol1.u ≈ sol2.u atol=1e-5

  # diagonal noise

  prob = SDEProblem(f,g,[u₀,u₀],tspan)
  sol = solve(prob,EulerHeun(),dt=dt,save_noise=true)

  _sol = deepcopy(sol) # to make sure the plot is correct
  W1 = NoiseGrid(reverse!(_sol.t),reverse!(_sol.W.W))
  prob1 = SDEProblem(f,g,sol[end],reverse(tspan),noise=W1)
  sol1 = solve(prob1,EulerHeun(),dt=dt)

  _sol = deepcopy(sol)
  W2 = NoiseWrapper(_sol.W, reverse=true)
  prob2 = SDEProblem(f,g,sol[end],reverse(tspan),noise=W2)
  sol2 = solve(prob2,EulerHeun(),dt=dt)

  @test sol.u ≈ reverse(sol1.u) atol=5e-2
  @test sol.u ≈ reverse(sol2.u) atol=5e-2
  @test sol1.u ≈ sol2.u atol=1e-5

end


@testset "Reverse a given NoiseProcess " begin
  # Noise Wrapper
  _W = WienerProcess(0.0,0.0,0.0)

  dt = 0.1
  calculate_step!(_W,dt,nothing,nothing)

  for i in 1:100
    accept_step!(_W,dt,nothing,nothing)
  end

  W1 = NoiseWrapper(_W, reverse=true)
  dt = -0.1
  calculate_step!(W1,dt,nothing,nothing)
  for i in 1:99
    accept_step!(W1,dt,nothing,nothing)
  end

  W2 = reverse(_W)
  dt = -0.1
  calculate_step!(W2,dt,nothing,nothing)
  for i in 1:99
    accept_step!(W2,dt,nothing,nothing)
  end

  @test isapprox(_W.W[2:end], reverse(W1.W), atol=1e-16)
  @test isapprox(_W.W[2:end], reverse(W2.W), atol=1e-16)

  # Noise Grid
  dt = 0.001
  t = 0:dt:1
  brownian_values = cumsum([0;[sqrt(dt)*randn() for i in 1:length(t)-1]])
  _W = NoiseGrid(t,brownian_values)

  prob = NoiseProblem(_W,(0.0,1.0))
  sol = solve(prob;dt=dt)

  W1 = NoiseGrid(reverse(sol.t),reverse(sol.W))
  prob = NoiseProblem(W1,(1.0,0.0))
  sol1 = solve(prob;dt=-dt)

  W2 = reverse(_W)
  prob = NoiseProblem(W2,(1.0,0.0))
  sol2 = solve(prob;dt=-dt)

  @test isapprox(sol.W, reverse(sol1.W), atol=1e-16)
  @test isapprox(sol.W, reverse(sol2.W), atol=1e-16)

end


@testset "SDE Ito Basic Reversal Tests" begin

  n = 100000
  T = 2.0
  dt = T/n
  x0 = 0.3

  b(u,p,t) =  sin(t) + cos(u)
  σ(u,p,t) = pi + atan(u)

  seed = 10
  Random.seed!(seed)
  W = [0.0; cumsum(sqrt(dt)*randn(n))]
  #using Plots; plot(W)

  p = nothing

  """
  For ito integrals
  """

  x = x0 # starting point
  xs = [x]
  t = 0.0
  ts = [0.0]
  for i in 1:n
    t, x
    #@show x, dt, (W[i+1] - W[i]), x+b(x,p,t)*dt, σ(x,p,t)*(W[i+1] - W[i])
    x += b(x,p,t)*dt + σ(x,p,t)*(W[i+1] - W[i]) # this is an ito integral
    t += dt
    push!(xs, x)
	push!(ts, t)

  end

  z = xs[end]

  #plt = plot(ts,xs)
  #plt = plot(xs)
  # Now reverse... t

  dσ_dx(u,p,t) = 1/(1 + u^2) # d(arctan(x))/dx
  z = xs[end]
  @show "starting point" z
  ys = [z]
  t = T
  cs = [0.0]
  for i in n:-1:1
    t, z
    cor =  1/2*dσ_dx(z,p,t)*σ(z,p,t)
    z -= (b(z,p,t)-0*2*cor)*dt + σ(z,p,t)*(W[i+1] - W[i]) # reverse ito integral
    t -= dt
    push!(ys, z)
    push!(cs, cor)
  end

  #plot!(ts,reverse(ys))
  # difference between forward and backward
  #plot(reverse(ys)-xs)
  # correction terms
  #plot(cs)
  #using DiffEqNoiseProcess, StochasticDiffEq


  W1 = NoiseGrid(ts,W)
  prob1 = SDEProblem(b,σ,x0,(0.0,2.0-1e-11),noise=W1)
  sol1 = solve(prob1,EM(false),dt=dt,adaptive=false)

  @test isapprox(xs, sol1.u, atol=1e-8)


  W1rev = NoiseGrid(reverse(ts),reverse(W))
  prob1 = SDEProblem(b,σ,sol1.u[end],(sol1.t[end],sol1.t[1]),noise=W1rev)
  sol2 = solve(prob1,EM(false),dt=dt,adaptive=false)

  # plot(ts,reverse(ys))
  # plot!(reverse(ts), sol2.u)
  @test isapprox(ys,sol2.u, atol=1e-6)
  @test !isapprox(sol1.u,reverse(sol2.u), atol=1e-0)

  bwrong(u,p,t) =  b(u,p,t) - 1//2*dσ_dx(u,p,t)*σ(u,p,t)
  W1rev = NoiseGrid(reverse(ts),reverse(W))
  prob1 = SDEProblem(bwrong,σ,sol1.u[end],(sol1.t[end],sol1.t[1]),noise=W1rev)
  sol2 = solve(prob1,EM(false),dt=dt,adaptive=false)

  @test !isapprox(ys,sol2.u, atol=1e-0)
  @test !isapprox(sol1.u,reverse(sol2.u), atol=1e-0)

  bcorrected(u,p,t) =  b(u,p,t) - 2*1//2*dσ_dx(u,p,t)*σ(u,p,t)
  W1rev = NoiseGrid(reverse(ts),reverse(W))
  prob1 = SDEProblem(bcorrected,σ,sol1.u[end],(sol1.t[end],sol1.t[1]),noise=W1rev)
  sol2 = solve(prob1,EM(false),dt=dt,adaptive=false)

  #plot(ts, sol1.u - reverse(sol2.u))
  @test !isapprox(ys,sol2.u, atol=1e-6)
  @test isapprox(sol1.u,reverse(sol2.u), atol=1e-0)

end



"""
Some more Ito reversals
"""


@testset "SDE Ito Reversal Tests" begin
  Random.seed!(100)
  α=1.0
  β=0.3

  dt = 1e-3
  tspan = (0.0,1.0)
  u₀=1/2

  tarray =  collect(tspan[1]:dt:tspan[2])

  f!(du,u,p,t) = du .= α*u
  g!(du,u,p,t) = du .= β*u


  prob = SDEProblem(f!,g!,[u₀],tspan)
  sol = solve(prob,EM(),dt=dt,save_noise=true, adaptive=false)

  fcorrected!(du,u,p,t) = du .= (α-β^2)*u

  _sol = deepcopy(sol) # to make sure the plot is correct
  W1 = NoiseGrid(reverse(_sol.t),reverse(_sol.W.W)
   # ,reverse(_sol.W.Z)
    )
  prob1 = SDEProblem(fcorrected!,g!,sol[end],reverse(tspan),noise=W1)
  sol1 = solve(prob1,EM(),dt=dt, adaptive=false)

  _sol = deepcopy(sol)
  W2 = NoiseWrapper(_sol.W, reverse=true)
  prob2 = SDEProblem(fcorrected!,g!,sol[end],reverse(tspan),noise=W2)
  sol2 = solve(prob2,EM(),dt=dt, save_noise=false, adaptive=false)

  # same time steps

  @test sol.u ≈ reverse(sol1.u) rtol=1e-2
  @test sol.u ≈ reverse(sol2.u) rtol=1e-2
  @test sol1.u ≈ sol2.u atol=1e-8


  f(u,p,t) = α*u
  g(u,p,t) = β*u

  prob = SDEProblem(f,g,u₀,tspan)
  sol = solve(prob,EM(),dt=dt,save_noise=true, adaptive=false)

  fcorrected(u,p,t) = (α-β^2)*u

  _sol = deepcopy(sol) # to make sure the plot is correct
  W1 = NoiseGrid(reverse(_sol.t),reverse(_sol.W.W)
    # ,reverse(_sol.W.Z)
   )
  prob1 = SDEProblem(fcorrected,g,sol[end],reverse(tspan),noise=W1)
  sol1 = solve(prob1,EM(),dt=dt, adaptive=false)

  _sol = deepcopy(sol)
  W2 = NoiseWrapper(_sol.W, reverse=true)
  prob2 = SDEProblem(fcorrected,g,sol[end],reverse(tspan),noise=W2)
  sol2 = solve(prob2,EM(),dt=dt, adaptive=false)

  @test sol.u ≈ reverse(sol1.u) rtol=1e-2
  @test sol.u ≈ reverse(sol2.u) rtol=1e-2
  @test sol1.u ≈ sol2.u atol=1e-7

  # using Plots; plt = plot(sol)
  # plot!(reverse(sol1.t),reverse(sol1[1,:]))
  # plot(vcat(sol.u - reverse(sol1.u) ...))
end
