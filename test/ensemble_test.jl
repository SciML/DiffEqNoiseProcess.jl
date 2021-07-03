using DiffEqNoiseProcess, SciMLBase
W = WienerProcess(0.0,1.0,1.0)
prob = NoiseProblem(W,(0.0,1.0))
ensemble_prob = EnsembleProblem(prob)
sol = solve(prob; dt=0.1)
sol = solve(ensemble_prob, EnsembleSerial(); trajectories=100, dt = 0.1)
