using StochasticDiffEq,  DiffEqBase, DiffEqNoiseProcess, Base.Test
f1 = (t,u) -> 1.01u
g1 = (t,u) -> 1.01u
dt = 1//2^(4)
prob1 = SDEProblem(f1,g1,1.0,(0.0,1.0))
sol1 = solve(prob1,EM(),dt=dt)

W2 = NoiseWrapper(sol1.W)

prob1 = SDEProblem(f1,g1,1.0,(0.0,1.0),noise=W2)

sol2 = solve(prob1,EM(),dt=dt)

@test sol1.u ≈ sol2.u

W3 = NoiseWrapper(sol1.W)
prob2 = SDEProblem(f1,g1,1.0,(0.0,1.0),noise=W3)

dt = 1//2^(5)
sol3 = solve(prob2,EM(),dt=dt)

using Plots
plot(sol1)
plot!(sol2)
plot!(sol3)

plot(sol1.W)
plot!(sol2.W)
plot!(sol3.W)
