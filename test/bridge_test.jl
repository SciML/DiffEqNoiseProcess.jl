@testset "Brownian Bridge" begin
    using DiffEqNoiseProcess, DiffEqBase, Test, Random, DiffEqBase.EnsembleAnalysis

    Random.seed!(100)
    W = BrownianBridge(0.0, 1.0, 0.0, 1.0, 0.0, 0.0)
    prob = NoiseProblem(W, (0.0, 1.0))
    ensemble_prob = EnsembleProblem(prob)
    @time sol = solve(ensemble_prob, dt = 0.1, trajectories = 100000)

    # Spot check the mean and the variance
    qs = 0:0.1:1
    for i in 2:10
        q = qs[i]
        @test ≈(timestep_mean(sol, i), q, atol = 1e-2)
        @test ≈(timestep_meanvar(sol, i)[2], (1 - q) * q, atol = 1e-2)
    end
    @test ≈(timestep_mean(sol, 1)[1], 0.0, atol = 1e-16)
    @test ≈(timestep_meanvar(sol, 1)[2], 0.0, atol = 1e-16)
    @test ≈(timestep_mean(sol, 11)[1], 1.0, atol = 1e-16)
    @test ≈(timestep_meanvar(sol, 11)[2], 0.0, atol = 1e-16)

    μ = 0.10500000000000001
    σ = 0.1
    t0 = 0.0
    tend = 10.0
    GB0 = 1.0
    GBend = 4.0

    W = GeometricBrownianBridge(μ, σ, t0, tend, GB0, GBend)
    prob = NoiseProblem(W, (t0, tend))
    ensemble_prob = EnsembleProblem(prob)
    @time sol = solve(ensemble_prob, dt = 1.0, trajectories = 100000)
    ts = t0:1.0:tend
    for i in 2:10
        t = ts[i]
        mean = log(GB0) + (t - t0) / (tend - t0) * (log(GBend) - log(GB0))
        var = (t - t0) / (tend - t0) * (tend - t) * σ^2

        m, v = timestep_meanvar(sol, i)

        @test ≈(m, exp(mean + var / 2), atol = 1e-2)
        @test ≈(v, (exp(var) - 1) * exp(2 * mean + var), atol = 1e-2)
    end

    Random.seed!(100)
    r = 100 # should be independent of the rate, so make it crazy
    rate(u, p, t) = r
    W = CompoundPoissonBridge(rate, 0.0, 1.0, 0.0, 1.0)
    prob = NoiseProblem(W, (0.0, 1.0))
    ensemble_prob = EnsembleProblem(prob)
    @time sol = solve(ensemble_prob, dt = 0.1, trajectories = 100000)

    # Spot check the mean and the variance
    qs = 0:0.1:1
    for i in 2:10
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
    for i in 2:8
        q = qs[i]
        @test ≈(timestep_mean(sol, i), q, atol = 1e-2)
        @test ≈(timestep_meanvar(sol, i)[2], (1 - q) * q, atol = 1e-2)
    end
    @test ≈(timestep_mean(sol, 1)[1], 0.0, atol = 1e-16)
    @test ≈(timestep_meanvar(sol, 1)[2], 0.0, atol = 1e-16)
    @test ≈(timestep_mean(sol, Int(2^(W.tree_depth) + 1))[1], W.W[end], atol = 1e-16)
    @test ≈(timestep_meanvar(sol, Int(2^(W.tree_depth) + 1))[2], 0.0, atol = 1e-16)
end

@testset "Scalar Ou-Bridge" begin
    using DiffEqNoiseProcess,
          DiffEqBase, Test, Random, DiffEqBase.EnsembleAnalysis, StatsBase, Statistics
    sols_forward = []
    dt = 0.125
    t_max = 20.0
    n = Int(t_max / dt + 2)
    st = []
    Random.seed!(1234)
    for i in 1:500_000
        local ou = OrnsteinUhlenbeckProcess(1.5, 0.25, 0.2, 0.0, 0.0)
        local prob = NoiseProblem(ou, (0.0, t_max))
        local s = solve(prob; dt = dt)
        push!(sols_forward, s.u)
        if i == 1
            st = s.t
        end
    end

    st2 = []
    sols_interpolated = []
    Random.seed!(1234)
    ou = OrnsteinUhlenbeckProcess(1.5, 0.25, 0.2, 0.0, 0.0)
    for i in 1:500_000
        local prob = NoiseProblem(ou, (0.0, t_max))
        local s2 = solve(prob; dt = t_max)
        s2.(st)
        push!(sols_interpolated, s2.u)
        if i == 1
            st2 = s2.t
        end
    end

    Random.seed!(1234)
    ou = OrnsteinUhlenbeckProcess(1.5, 0.25, 0.2, 0.0, 0.0)
    prob = NoiseProblem(ou, (0.0, t_max))
    sols_solver = []
    st = []
    for i in 1:500_000
        local s = solve(prob; dt = t_max)
        ou_bridge = OrnsteinUhlenbeckBridge(1.5, 0.25, 0.2, 0.0, t_max, 0.0, s.u[end])
        s = solve(NoiseProblem(ou_bridge, (0.0, t_max)); dt = 0.125)
        push!(sols_solver, s.u)
        if i == 1
            st = s.t
        end
    end

    @test all(isapprox.(mean(sols_forward), mean(sols_interpolated); atol = 1e-3))
    @test all(isapprox.(mean(sols_forward), mean(sols_solver); atol = 1e-3))

    @test all(isapprox.(std(sols_forward), std(sols_interpolated); atol = 1e-3))
    @test all(isapprox.(std(sols_forward), std(sols_solver); atol = 1e-3))
end

@testset "Vector Ou-Bridge" begin
    sols_forward = []
    dt = 0.125
    t_max = 20.0
    n = Int(t_max / dt + 2)
    st = []
    Random.seed!(1234)
    ou = OrnsteinUhlenbeckProcess([1.5], [0.25], [0.2], 0.0, [0.0])
    for i in 1:500_000
        local prob = NoiseProblem(ou, (0.0, t_max))
        local s = solve(prob; dt = dt)
        push!(sols_forward, s.u)
        if i == 1
            st = s.t
        end
    end

    st2 = []
    sols_interpolated = []
    Random.seed!(1234)
    ou = OrnsteinUhlenbeckProcess([1.5], [0.25], [0.2], 0.0, [0.0])
    for i in 1:500_000
        local prob = NoiseProblem(ou, (0.0, t_max))
        local s2 = solve(prob; dt = t_max)
        s2.(st)
        push!(sols_interpolated, s2.u)
        if i == 1
            st2 = s2.t
        end
    end

    Random.seed!(1234)
    ou = OrnsteinUhlenbeckProcess([1.5], [0.25], [0.2], 0.0, [0.0])
    prob = NoiseProblem(ou, (0.0, t_max))
    sols_solver = []
    st = []
    for i in 1:500_000
        local s = solve(prob; dt = t_max)
        local ou_bridge = OrnsteinUhlenbeckBridge([1.5],
            [0.25],
            [0.2],
            0.0,
            t_max,
            [0.0],
            s.u[end])
        local s = solve(NoiseProblem(ou_bridge, (0.0, t_max)); dt = 0.125)
        push!(sols_solver, s.u)
        if i == 1
            st = s.t
        end
    end

    sols_forward = map(sols_forward) do x
        [y[1] for y in x]
    end
    sols_interpolated = map(sols_interpolated) do x
        [y[1] for y in x]
    end
    sols_solver = map(sols_solver) do x
        [y[1] for y in x]
    end

    @test all(isapprox.(mean(sols_forward), mean(sols_interpolated); atol = 1e-3))
    @test all(isapprox.(mean(sols_forward), mean(sols_solver); atol = 1e-3))

    @test all(isapprox.(std(sols_forward), std(sols_interpolated); atol = 1e-3))
    @test all(isapprox.(std(sols_forward), std(sols_solver); atol = 1e-3))
end

@testset "Inplace OU-Bridge" begin
    sols_forward = []
    dt = 0.125
    t_max = 20.0
    n = Int(t_max / dt + 2)
    st = []
    Random.seed!(1234)
    ou = OrnsteinUhlenbeckProcess!([1.5], [0.25], [0.2], 0.0, [0.0])
    for i in 1:500_000
        local prob = NoiseProblem(ou, (0.0, t_max))
        local s = solve(prob; dt = dt)
        push!(sols_forward, s.u)
        if i == 1
            st = s.t
        end
    end

    st2 = []
    sols_interpolated = []
    Random.seed!(1234)
    ou = OrnsteinUhlenbeckProcess!([1.5], [0.25], [0.2], 0.0, [0.0])
    for i in 1:500_000
        local prob = NoiseProblem(ou, (0.0, t_max))
        local s2 = solve(prob; dt = t_max)
        s2.(st)
        push!(sols_interpolated, s2.u)
        if i == 1
            st2 = s2.t
        end
    end

    Random.seed!(1234)
    ou = OrnsteinUhlenbeckProcess!([1.5], [0.25], [0.2], 0.0, [0.0])
    prob = NoiseProblem(ou, (0.0, t_max))
    sols_solver = []
    st = []
    for i in 1:500_000
        local s = solve(prob; dt = t_max)
        local ou_bridge = OrnsteinUhlenbeckBridge!([1.5],
            [0.25],
            [0.2],
            0.0,
            t_max,
            [0.0],
            s.u[end])
        local s = solve(NoiseProblem(ou_bridge, (0.0, t_max)); dt = 0.125)
        push!(sols_solver, s.u)
        if i == 1
            st = s.t
        end
    end

    sols_forward = map(sols_forward) do x
        [y[1] for y in x]
    end
    sols_interpolated = map(sols_interpolated) do x
        [y[1] for y in x]
    end
    sols_solver = map(sols_solver) do x
        [y[1] for y in x]
    end

    @test all(isapprox.(mean(sols_forward), mean(sols_interpolated); atol = 1e-3))
    @test all(isapprox.(mean(sols_forward), mean(sols_solver); atol = 1e-3))

    @test all(isapprox.(std(sols_forward), std(sols_interpolated); atol = 1e-3))
    @test all(isapprox.(std(sols_forward), std(sols_solver); atol = 1e-3))
end
