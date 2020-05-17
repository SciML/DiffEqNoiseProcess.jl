using StochasticDiffEq, DiffEqNoiseProcess, Test
α=1
β=1
u₀=1/2
f(u,p,t) = α*u
g(u,p,t) = β*u
dt = 1//2^(4)
tspan = (0.0,1.0)
prob = SDEProblem(f,g,u₀,(0.0,1.0))
sol =solve(prob,EulerHeun(),dt=0.01,save_noise=true)
_sol = deepcopy(sol) # to make sure the plot is correct
W3 = NoiseGrid(reverse!(_sol.t),reverse!(_sol.W))
prob3 = SDEProblem(f,g,sol[end],(1.0,0.0),noise=W3)
sol2 = solve(prob3,EulerHeun(),dt=0.01)
@test sol.u ≈ reverse!(sol2.u) atol=1e-1
