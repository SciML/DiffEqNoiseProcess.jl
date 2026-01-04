@testset "Reinit" begin
    using StochasticDiffEq, DiffEqNoiseProcess, Statistics

    f = (u, p, t) -> 1.0
    g = (u, p, t) -> 0.2
    tspan = (0.5, 1.0)

    @testset "WienerProcess" begin
        W = WienerProcess(0.0, 0.0)
        t = 1.0
        expected_mean = 0.0
        expected_variance = t^2

        meanW = mean(W.curW for i in 1:10_000 if reinit!(W, 0.0, t0 = t) === nothing)
        varW = var(W.curW for i in 1:10_000 if reinit!(W, 0.0, t0 = t) === nothing)
        @test meanW ≈ expected_mean atol = 0.1
        @test varW ≈ expected_variance rtol = 0.1

        prob = NoiseProblem(W, tspan)
        ensemble_prob = EnsembleProblem(prob, output_func = (sol, i) -> (sol.u[end], false))
        sol = solve(ensemble_prob, dt = 1 / 10, trajectories = 40_000)
        @test mean(sol) ≈ expected_mean atol = 0.1
        @test var(sol) ≈ expected_variance atol = 0.1

        prob = SDEProblem(f, g, 1.0, tspan, noise = W, save_noise = true)

        ensemble_probW = EnsembleProblem(
            prob,
            output_func = (sol, i) -> (sol.W.W[end], false)
        )

        solW_at_1 = solve(ensemble_probW, EM(), dt = 1 / 10, trajectories = 40_000)

        @test mean(solW_at_1) ≈ expected_mean atol = 0.1
        @test var(solW_at_1) ≈ expected_variance atol = 0.1
    end

    @testset "GBM" begin
        μ = 0.5
        σ = 0.1
        W0 = 1.0
        W = GeometricBrownianMotionProcess(μ, σ, 0.0, W0, nothing)

        t = 1.0
        expected_mean = W0 * exp(μ * t)
        expected_variance = W0^2 * exp(2μ * t) * (exp(σ^2 * t) - 1)

        meanW = mean(W.curW for i in 1:10_000 if reinit!(W, 0.0, t0 = t) === nothing)
        varW = var(W.curW for i in 1:10_000 if reinit!(W, 0.0, t0 = t) === nothing)
        @test meanW ≈ expected_mean atol = 0.1
        @test varW ≈ expected_variance rtol = 0.1

        prob = NoiseProblem(W, tspan)
        ensemble_prob = EnsembleProblem(prob, output_func = (sol, i) -> (sol.u[end], false))
        sol = solve(ensemble_prob, dt = 1 / 10, trajectories = 40_000)
        @test mean(sol) ≈ expected_mean rtol = 0.1
        @test var(sol) ≈ expected_variance atol = 0.1

        prob = SDEProblem(f, g, 1.0, tspan, noise = W, save_noise = true)

        ensemble_probW = EnsembleProblem(
            prob,
            output_func = (sol, i) -> (sol.W.W[end], false)
        )

        solW_at_1 = solve(ensemble_probW, EM(), dt = 1 / 10, trajectories = 40_000)

        @test mean(solW_at_1) ≈ expected_mean rtol = 0.1
        @test var(solW_at_1) ≈ expected_variance atol = 0.1
    end

    @testset "OrnsteinUhlenbeck" begin
        Θ = 1.0
        μ = 1.2
        σ = 0.3
        W0 = 1.0
        W = OrnsteinUhlenbeckProcess(Θ, μ, σ, 0.0, W0, nothing)

        t = 1.0
        expected_mean = μ + (W0 - μ) * exp(-Θ * t)
        expected_variance = (1 - exp(-2Θ * t)) * σ^2 / (2Θ)

        meanW = mean(W.curW for i in 1:10_000 if reinit!(W, 0.0, t0 = t) === nothing)
        varW = var(W.curW for i in 1:10_000 if reinit!(W, 0.0, t0 = t) === nothing)
        @test meanW ≈ expected_mean atol = 0.1
        @test varW ≈ expected_variance rtol = 0.1

        prob = NoiseProblem(W, tspan)
        ensemble_prob = EnsembleProblem(prob, output_func = (sol, i) -> (sol.u[end], false))
        sol = solve(ensemble_prob, dt = 1 / 10, trajectories = 40_000)
        @test mean(sol) ≈ expected_mean rtol = 0.1
        @test var(sol) ≈ expected_variance atol = 0.1

        prob = SDEProblem(f, g, 1.0, tspan, noise = W, save_noise = true)

        ensemble_probW = EnsembleProblem(
            prob,
            output_func = (sol, i) -> (sol.W.W[end], false)
        )

        solW_at_1 = solve(ensemble_probW, EM(), dt = 1 / 10, trajectories = 40_000)

        @test mean(solW_at_1) ≈ expected_mean rtol = 0.1
        @test var(solW_at_1) ≈ expected_variance atol = 0.1
    end

    @testset "BrownianBridge" begin
        W = BrownianBridge(0.0, 1.0, 0.0, 1.0, nothing, nothing)

        t = 1.0
        expected_mean = 1.0
        expected_variance = 0.0

        meanW = mean(W.curW for i in 1:10_000 if reinit!(W, 0.0, t0 = t) === nothing)
        varW = var(W.curW for i in 1:10_000 if reinit!(W, 0.0, t0 = t) === nothing)
        @test meanW ≈ expected_mean rtol = 0.1
        @test varW ≈ expected_variance atol = 0.1

        prob = NoiseProblem(W, tspan)
        ensemble_prob = EnsembleProblem(prob, output_func = (sol, i) -> (sol.u[end], false))
        sol = solve(ensemble_prob, dt = 1 / 10, trajectories = 40_000)
        @test mean(sol) ≈ expected_mean rtol = 0.1
        @test var(sol) ≈ expected_variance atol = 0.1

        prob = SDEProblem(f, g, 1.0, tspan, noise = W, save_noise = true)

        ensemble_probW = EnsembleProblem(
            prob,
            output_func = (sol, i) -> (sol.W.W[end], false)
        )

        solW_at_1 = solve(ensemble_probW, EM(), dt = 1 / 10, trajectories = 40_000)

        @test mean(solW_at_1) ≈ expected_mean rtol = 0.1
        @test var(solW_at_1) ≈ expected_variance atol = 0.1
    end

    @testset "GBMBridge" begin
        μ, σ = 1.2, 0.2
        t0, tend = 0.0, 1.0
        W0, Wend = 1.0, 2.0
        Z0 = Zend = nothing
        W = GeometricBrownianBridge(μ, σ, t0, tend, W0, Wend, Z0, Zend)

        t = 1.0
        expected_mean = Wend
        expected_variance = 0.0

        meanW = mean(W.curW for i in 1:10_000 if reinit!(W, 0.0, t0 = t) === nothing)
        varW = var(W.curW for i in 1:10_000 if reinit!(W, 0.0, t0 = t) === nothing)
        @test meanW ≈ expected_mean rtol = 0.1
        @test varW ≈ expected_variance atol = 0.1

        prob = NoiseProblem(W, tspan)
        ensemble_prob = EnsembleProblem(prob, output_func = (sol, i) -> (sol.u[end], false))
        sol = solve(ensemble_prob, dt = 1 / 10, trajectories = 40_000)
        @test mean(sol) ≈ expected_mean rtol = 0.1
        @test var(sol) ≈ expected_variance atol = 0.1

        prob = SDEProblem(f, g, 1.0, tspan, noise = W, save_noise = true)

        ensemble_probW = EnsembleProblem(
            prob,
            output_func = (sol, i) -> (sol.W.W[end], false)
        )

        solW_at_1 = solve(ensemble_probW, EM(), dt = 1 / 10, trajectories = 40_000)

        @test mean(solW_at_1) ≈ expected_mean rtol = 0.1
        @test var(solW_at_1) ≈ expected_variance atol = 0.1
    end

    @testset "CompoundPoissonBridge" begin
        rate = (u, p, t) -> 100
        W = CompoundPoissonBridge(rate, 0.0, 1.0, 0.0, 1.0)

        t = 1.0
        expected_mean = 1.0
        expected_variance = 0.0

        meanW = mean(W.curW for i in 1:10_000 if reinit!(W, 0.0, t0 = t) === nothing)
        varW = var(W.curW for i in 1:10_000 if reinit!(W, 0.0, t0 = t) === nothing)
        @test meanW ≈ expected_mean rtol = 0.1
        @test varW ≈ expected_variance atol = 0.1

        prob = NoiseProblem(W, tspan)
        ensemble_prob = EnsembleProblem(prob, output_func = (sol, i) -> (sol.u[end], false))
        sol = solve(ensemble_prob, dt = 1 / 10, trajectories = 40_000)
        @test mean(sol) ≈ expected_mean rtol = 0.1
        @test var(sol) ≈ expected_variance atol = 0.1

        prob = SDEProblem(f, g, 1.0, tspan, noise = W, save_noise = true)

        ensemble_probW = EnsembleProblem(
            prob,
            output_func = (sol, i) -> (sol.W.W[end], false)
        )

        solW_at_1 = solve(ensemble_probW, EM(), dt = 1 / 10, trajectories = 40_000)

        @test mean(solW_at_1) ≈ expected_mean rtol = 0.1
        @test var(solW_at_1) ≈ expected_variance atol = 0.1
    end

    @testset "NoiseTransport" begin
        W = NoiseTransport(0.0, (u, p, t, Y) -> Y * t, randn)
        t = 2.0
        expected_mean = 0.0
        expected_variance = t^2

        meanW = mean(W.curW for i in 1:10_000 if reinit!(W, 0.0, t0 = t) === nothing)
        varW = var(W.curW for i in 1:10_000 if reinit!(W, 0.0, t0 = t) === nothing)
        @test meanW ≈ expected_mean atol = 0.1
        @test varW ≈ expected_variance rtol = 0.1

        prob = NoiseProblem(W, tspan)
        ensemble_prob = EnsembleProblem(
            prob,
            output_func = (sol, i) -> (first(sol(t)), false)
        )
        sol = solve(ensemble_prob, dt = 1 / 10, trajectories = 40_000)

        @test mean(sol) ≈ expected_mean atol = 0.1
        @test var(sol) ≈ expected_variance atol = 0.1

        prob = SDEProblem(f, g, 1.0, tspan, noise = W, save_noise = true)

        ensemble_probW = EnsembleProblem(
            prob,
            output_func = (sol, i) -> (first(sol.W(t)), false)
        )
        solW_at_1 = solve(ensemble_probW, EM(), dt = 1 / 10, trajectories = 40_000)

        @test mean(solW_at_1) ≈ expected_mean atol = 0.1
        @test var(solW_at_1) ≈ expected_variance atol = 0.1
    end
end
