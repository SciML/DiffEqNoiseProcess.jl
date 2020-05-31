@testset "SDE Reversal Tests" begin

using StochasticDiffEq, DiffEqNoiseProcess, Test, Random
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
sol =solve(prob,EulerHeun(),dt=dt,save_noise=true, adaptive=false)

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
