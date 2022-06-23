@testset "Brownian Bridge" begin

    using DiffEqNoiseProcess, DiffEqBase, Test, Random, DiffEqBase.EnsembleAnalysis

    Random.seed!(100)
    W = BrownianBridge(0.0, 1.0, 0.0, 1.0, 0.0, 0.0)
    prob = NoiseProblem(W, (0.0, 1.0))
    ensemble_prob = EnsembleProblem(prob)
    @time sol = solve(ensemble_prob, dt = 0.1, trajectories = 100000)

    # Spot check the mean and the variance
    qs = 0:0.1:1
    for i = 2:10
        q = qs[i]
        @test ≈(timestep_mean(sol, i), q, atol = 1e-2)
        @test ≈(timestep_meanvar(sol, i)[2], (1 - q) * q, atol = 1e-2)
    end
    @test ≈(timestep_mean(sol, 1)[1], 0.0, atol = 1e-16)
    @test ≈(timestep_meanvar(sol, 1)[2], 0.0, atol = 1e-16)
    @test ≈(timestep_mean(sol, 11)[1], 1.0, atol = 1e-16)
    @test ≈(timestep_meanvar(sol, 11)[2], 0.0, atol = 1e-16)


    μ = 1.2
    σ = 2.2
    W = GeometricBrownianBridge(μ, σ, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0)
    prob = NoiseProblem(W, (0.0, 1.0))
    ensemble_prob = EnsembleProblem(prob)
    @time sol = solve(ensemble_prob, dt = 0.1, trajectories = 100)

    Random.seed!(100)
    r = 100 # should be independent of the rate, so make it crazy
    rate(u, p, t) = r
    W = CompoundPoissonBridge(rate, 0.0, 1.0, 0.0, 1.0)
    prob = NoiseProblem(W, (0.0, 1.0))
    ensemble_prob = EnsembleProblem(prob)
    @time sol = solve(ensemble_prob, dt = 0.1, trajectories = 100000)

    # Spot check the mean and the variance
    qs = 0:0.1:1
    for i = 2:10
        q = qs[i]
        # Mean and variance of binomial matches that of the Brownian bridge!
        @test ≈(timestep_mean(sol, i), q, atol = 1e-2)
        @test ≈(timestep_meanvar(sol, i)[2], (1 - q) * q, atol = 1e-2)
    end
    @test ≈(timestep_mean(sol, 1)[1], 0.0, atol = 1e-16)
    @test ≈(timestep_meanvar(sol, 1)[2], 0.0, atol = 1e-16)
    @test ≈(timestep_mean(sol, 11)[1], 1.0, atol = 1e-16)
    @test ≈(timestep_meanvar(sol, 11)[2], 0.0, atol = 1e-16)

    # check VBT distributional properties

    W = VirtualBrownianTree(0.0, 0.0; Wend = 1.0, tree_depth = 3)
    prob = NoiseProblem(W, (0.0, 1.0))
    function prob_func(prob, i, repeat)
        # to re-instantiate PRNG
        Wtmp = VirtualBrownianTree(0.0, 0.0; Wend = 1.0, tree_depth = 3)
        remake(prob, noise = Wtmp)
    end
    ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
    @time sol = solve(ensemble_prob, dt = 0.125, trajectories = 100000)

    # Spot check the mean and the variance
    qs = 0:0.125:1
    for i = 2:8
        q = qs[i]
        @test ≈(timestep_mean(sol, i), q, atol = 1e-2)
        @test ≈(timestep_meanvar(sol, i)[2], (1 - q) * q, atol = 1e-2)
    end
    @test ≈(timestep_mean(sol, 1)[1], 0.0, atol = 1e-16)
    @test ≈(timestep_meanvar(sol, 1)[2], 0.0, atol = 1e-16)
    @test ≈(timestep_mean(sol, Int(2^(W.tree_depth) + 1))[1], W.W[end], atol = 1e-16)
    @test ≈(timestep_meanvar(sol, Int(2^(W.tree_depth) + 1))[2], 0.0, atol = 1e-16)

end
