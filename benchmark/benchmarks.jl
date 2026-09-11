using DiffEqNoiseProcess, DiffEqBase, BenchmarkTools

const SUITE = BenchmarkGroup()

# =============================================================================
# Noise process simulation
# =============================================================================

SUITE["simulate"] = BenchmarkGroup()

W_scalar = WienerProcess(0.0, 0.0, 0.0)
prob_w = NoiseProblem(W_scalar, (0.0, 1.0); seed = 1234)

W_vec = WienerProcess!(0.0, zeros(10), zeros(10))
prob_wv = NoiseProblem(W_vec, (0.0, 1.0); seed = 1234)

gbm = GeometricBrownianMotionProcess(0.01, 0.2, 0.0, 1.0, 1.0)
prob_gbm = NoiseProblem(gbm, (0.0, 1.0); seed = 1234)

ou = OrnsteinUhlenbeckProcess(1.0, 0.0, 0.5, 0.0, 0.0)
prob_ou = NoiseProblem(ou, (0.0, 1.0); seed = 1234)

SUITE["simulate"]["wiener_scalar"] = @benchmarkable solve($prob_w; dt = 0.001)
SUITE["simulate"]["wiener_vector"] = @benchmarkable solve($prob_wv; dt = 0.001)
SUITE["simulate"]["geometric_brownian"] = @benchmarkable solve($prob_gbm; dt = 0.001)
SUITE["simulate"]["ornstein_uhlenbeck"] = @benchmarkable solve($prob_ou; dt = 0.001)

# =============================================================================
# Bridge processes
# =============================================================================

SUITE["bridge"] = BenchmarkGroup()

Wb = BrownianBridge(0.0, 1.0, 0.0, 1.0)
prob_wb = NoiseProblem(Wb, (0.0, 1.0); seed = 1234)

gbmb = GeometricBrownianBridge(0.01, 0.2, 0.0, 1.0, 1.0, 1.5)
prob_gbmb = NoiseProblem(gbmb, (0.0, 1.0); seed = 1234)

SUITE["bridge"]["brownian"] = @benchmarkable solve($prob_wb; dt = 0.001)
SUITE["bridge"]["geometric_brownian"] = @benchmarkable solve($prob_gbmb; dt = 0.001)

# =============================================================================
# Ensemble simulation
# =============================================================================

SUITE["ensemble"] = BenchmarkGroup()

eprob = EnsembleProblem(prob_w)
SUITE["ensemble"]["trajectories_100"] = @benchmarkable solve(
    $eprob; dt = 0.01, trajectories = 100
)
