using DiffEqNoiseProcess, DiffEqBase, DiffEqMonteCarlo,
      Base.Test, DataStructures

W = BrownianBridge(0.0,1.0,0.0,0.0,0.0,0.0)
dt = 0.1
W.dt = dt
DiffEqNoiseProcess.setup_next_step!(W)
for i in 1:10
  accept_step!(W,dt)
end
@test ≈(W[end],0.0,atol=3e-8)

W = BrownianBridge(0.0,1.0,0.0,0.0,0.0,0.0)
prob = NoiseProblem(W,(0.0,1.0))
sol = solve(prob;dt=0.1)

@test ≈(sol[end],0.0,atol=3e-8)

W = BrownianBridge(0.0,1.0,0.0,1.0,0.0,0.0)
prob = NoiseProblem(W,(0.0,1.0))
sol = solve(prob;dt=0.1)

@test ≈(sol[end],1.0,atol=3e-8)

#=
const μ = 1.0
const σ = 2.0

W = BrownianBridge(0.0,1.0,0.0,0.0,0.0,0.0,rswm=RSWM(adaptivealg=:RSwM1))
prob = NoiseProblem(W,(0.0,1.0))
sol = solve(prob;dt=0.1)

using Plots; plot(sol)

sol.t
sol[end]

W.S₁
=#
