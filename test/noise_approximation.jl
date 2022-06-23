@testset "NoiseApproximation" begin

    using DiffEqNoiseProcess, DiffEqBase, StochasticDiffEq
    using Test

    using DiffEqProblemLibrary.SDEProblemLibrary: importsdeproblems
    importsdeproblems()
    import DiffEqProblemLibrary.SDEProblemLibrary: prob_sde_linear, prob_sde_2Dlinear

    prob = prob_sde_linear
    integrator = init(prob, EM(), dt = 0.01)

    W = NoiseApproximation(integrator)


    dt = 0.1
    calculate_step!(W, dt, nothing, nothing)
    dWold = W.dW
    @test W.curW == W[1]
    @test W.curt == 0.0
    accept_step!(W, dt, nothing, nothing)
    @test W.curW == 0.5 + dWold
    @test W.curt == dt
    @test W.curW + W.dW == W[end]
    @test W.curW == W[11]

    W = NoiseApproximation(integrator)
    for i = 1:10
        accept_step!(W, dt, nothing, nothing)
    end
    W.t[end] == 1.0

    W = NoiseApproximation(integrator)
    prob = NoiseProblem(W, (0.0, 1.0))
    sol = solve(prob; dt = 0.1)


    prob = prob_sde_2Dlinear
    integrator = init(prob, EM(), dt = 0.01)
    W = NoiseApproximation(integrator)
    prob = NoiseProblem(W, (0.0, 1.0))
    sol = solve(prob; dt = 0.1)

end
