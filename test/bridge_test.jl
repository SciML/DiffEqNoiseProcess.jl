# Memory discipline for this file: it draws 4 x 100k ensemble trajectories and
# 9 x 500k bridge trajectories. Retaining full solution objects peaks at ~21GB
# resident and OOMs the 7GB hosted CI runners, so the ensembles keep only the
# value arrays (timestep_mean/timestep_meanvar need `el.u[i]`) and the Ou-Bridge
# testsets accumulate streaming sums instead of solution collections — only the
# elementwise mean/std are asserted.
keep_u = (sol, ctx) -> ((u = sol.u,), false)

mutable struct RunningStats
    n::Int
    sum::Vector{Float64}
    sumsq::Vector{Float64}
    RunningStats() = new(0, Float64[], Float64[])
end
function push_stats!(rs::RunningStats, u)
    if rs.n == 0
        resize!(rs.sum, length(u))
        fill!(rs.sum, 0.0)
        resize!(rs.sumsq, length(u))
        fill!(rs.sumsq, 0.0)
    end
    rs.n += 1
    for i in eachindex(u)
        v = first(u[i])
        rs.sum[i] += v
        rs.sumsq[i] += abs2(v)
    end
    return rs
end
stats_mean(rs::RunningStats) = rs.sum ./ rs.n
function stats_std(rs::RunningStats)
    m = stats_mean(rs)
    return sqrt.(max.(rs.sumsq ./ rs.n .- abs2.(m), 0.0) .* (rs.n / (rs.n - 1)))
end

@testset "Brownian Bridge" begin
    using DiffEqNoiseProcess, DiffEqBase, Test, Random, DiffEqBase.EnsembleAnalysis

    Random.seed!(100)
    W = BrownianBridge(0.0, 1.0, 0.0, 1.0, 0.0, 0.0)
    prob = NoiseProblem(W, (0.0, 1.0))
    ensemble_prob = EnsembleProblem(prob, output_func = keep_u)
    @time sol = solve(ensemble_prob, dt = 0.1, trajectories = 100000)

    # Spot check the mean and the variance
    qs = 0:0.1:1
    for i in 2:10
        q = qs[i]
        @test ≈(timestep_mean(sol, i), q, atol = 1.0e-2)
        @test ≈(timestep_meanvar(sol, i)[2], (1 - q) * q, atol = 1.0e-2)
    end
    @test ≈(timestep_mean(sol, 1)[1], 0.0, atol = 1.0e-16)
    @test ≈(timestep_meanvar(sol, 1)[2], 0.0, atol = 1.0e-16)
    @test ≈(timestep_mean(sol, 11)[1], 1.0, atol = 1.0e-16)
    @test ≈(timestep_meanvar(sol, 11)[2], 0.0, atol = 1.0e-16)

    μ = 0.10500000000000001
    σ = 0.1
    t0 = 0.0
    tend = 10.0
    GB0 = 1.0
    GBend = 4.0

    W = GeometricBrownianBridge(μ, σ, t0, tend, GB0, GBend)
    prob = NoiseProblem(W, (t0, tend))
    ensemble_prob = EnsembleProblem(prob, output_func = keep_u)
    @time sol = solve(ensemble_prob, dt = 1.0, trajectories = 100000)
    ts = t0:1.0:tend
    for i in 2:10
        t = ts[i]
        mean = log(GB0) + (t - t0) / (tend - t0) * (log(GBend) - log(GB0))
        var = (t - t0) / (tend - t0) * (tend - t) * σ^2

        m, v = timestep_meanvar(sol, i)

        @test ≈(m, exp(mean + var / 2), atol = 1.0e-2)
        @test ≈(v, (exp(var) - 1) * exp(2 * mean + var), atol = 1.0e-2)
    end

    Random.seed!(100)
    r = 100 # should be independent of the rate, so make it crazy
    rate(u, p, t) = r
    W = CompoundPoissonBridge(rate, 0.0, 1.0, 0.0, 1.0)
    prob = NoiseProblem(W, (0.0, 1.0))
    ensemble_prob = EnsembleProblem(prob, output_func = keep_u)
    @time sol = solve(ensemble_prob, dt = 0.1, trajectories = 100000)

    # Spot check the mean and the variance
    qs = 0:0.1:1
    for i in 2:10
        q = qs[i]
        # Mean and variance of binomial matches that of the Brownian bridge!
        @test ≈(timestep_mean(sol, i), q, atol = 1.0e-2)
        @test ≈(timestep_meanvar(sol, i)[2], (1 - q) * q, atol = 1.0e-2)
    end
    @test ≈(timestep_mean(sol, 1)[1], 0.0, atol = 1.0e-16)
    @test ≈(timestep_meanvar(sol, 1)[2], 0.0, atol = 1.0e-16)
    @test ≈(timestep_mean(sol, 11)[1], 1.0, atol = 1.0e-16)
    @test ≈(timestep_meanvar(sol, 11)[2], 0.0, atol = 1.0e-16)

    # check VBT distributional properties

    W = VirtualBrownianTree(0.0, 0.0; Wend = 1.0, tree_depth = 3)
    prob = NoiseProblem(W, (0.0, 1.0))
    function prob_func(prob, ctx)
        # to re-instantiate PRNG
        Wtmp = VirtualBrownianTree(0.0, 0.0; Wend = 1.0, tree_depth = 3)
        remake(prob, noise = Wtmp)
    end
    ensemble_prob = EnsembleProblem(prob, prob_func = prob_func, output_func = keep_u)
    @time sol = solve(ensemble_prob, dt = 0.125, trajectories = 100000)

    # Spot check the mean and the variance
    qs = 0:0.125:1
    for i in 2:8
        q = qs[i]
        @test ≈(timestep_mean(sol, i), q, atol = 1.0e-2)
        @test ≈(timestep_meanvar(sol, i)[2], (1 - q) * q, atol = 1.0e-2)
    end
    @test ≈(timestep_mean(sol, 1)[1], 0.0, atol = 1.0e-16)
    @test ≈(timestep_meanvar(sol, 1)[2], 0.0, atol = 1.0e-16)
    @test ≈(timestep_mean(sol, Int(2^(W.tree_depth) + 1))[1], W.W[end], atol = 1.0e-16)
    @test ≈(timestep_meanvar(sol, Int(2^(W.tree_depth) + 1))[2], 0.0, atol = 1.0e-16)
end

@testset "Scalar Ou-Bridge" begin
    using DiffEqNoiseProcess,
        DiffEqBase, Test, Random, DiffEqBase.EnsembleAnalysis, StatsBase, Statistics
    stats_forward = RunningStats()
    dt = 0.125
    t_max = 20.0
    st = []
    Random.seed!(1234)
    for i in 1:500_000
        local ou = OrnsteinUhlenbeckProcess(1.5, 0.25, 0.2, 0.0, 0.0)
        local prob = NoiseProblem(ou, (0.0, t_max))
        local s = solve(prob; dt = dt)
        push_stats!(stats_forward, s.u)
        if i == 1
            st = s.t
        end
    end

    stats_interpolated = RunningStats()
    Random.seed!(1234)
    ou = OrnsteinUhlenbeckProcess(1.5, 0.25, 0.2, 0.0, 0.0)
    for i in 1:500_000
        local prob = NoiseProblem(ou, (0.0, t_max))
        local s2 = solve(prob; dt = t_max)
        s2.(st)
        push_stats!(stats_interpolated, s2.u)
    end

    Random.seed!(1234)
    ou = OrnsteinUhlenbeckProcess(1.5, 0.25, 0.2, 0.0, 0.0)
    prob = NoiseProblem(ou, (0.0, t_max))
    stats_solver = RunningStats()
    for i in 1:500_000
        local s = solve(prob; dt = t_max)
        ou_bridge = OrnsteinUhlenbeckBridge(1.5, 0.25, 0.2, 0.0, t_max, 0.0, s.u[end])
        s = solve(NoiseProblem(ou_bridge, (0.0, t_max)); dt = 0.125)
        push_stats!(stats_solver, s.u)
    end

    @test all(isapprox.(stats_mean(stats_forward), stats_mean(stats_interpolated); atol = 1.0e-3))
    @test all(isapprox.(stats_mean(stats_forward), stats_mean(stats_solver); atol = 1.0e-3))

    @test all(isapprox.(stats_std(stats_forward), stats_std(stats_interpolated); atol = 1.0e-3))
    @test all(isapprox.(stats_std(stats_forward), stats_std(stats_solver); atol = 1.0e-3))
end

@testset "Vector Ou-Bridge" begin
    stats_forward = RunningStats()
    dt = 0.125
    t_max = 20.0
    st = []
    Random.seed!(1234)
    ou = OrnsteinUhlenbeckProcess([1.5], [0.25], [0.2], 0.0, [0.0])
    for i in 1:500_000
        local prob = NoiseProblem(ou, (0.0, t_max))
        local s = solve(prob; dt = dt)
        push_stats!(stats_forward, s.u)
        if i == 1
            st = s.t
        end
    end

    stats_interpolated = RunningStats()
    Random.seed!(1234)
    ou = OrnsteinUhlenbeckProcess([1.5], [0.25], [0.2], 0.0, [0.0])
    for i in 1:500_000
        local prob = NoiseProblem(ou, (0.0, t_max))
        local s2 = solve(prob; dt = t_max)
        s2.(st)
        push_stats!(stats_interpolated, s2.u)
    end

    Random.seed!(1234)
    ou = OrnsteinUhlenbeckProcess([1.5], [0.25], [0.2], 0.0, [0.0])
    prob = NoiseProblem(ou, (0.0, t_max))
    stats_solver = RunningStats()
    for i in 1:500_000
        local s = solve(prob; dt = t_max)
        local ou_bridge = OrnsteinUhlenbeckBridge(
            [1.5],
            [0.25],
            [0.2],
            0.0,
            t_max,
            [0.0],
            s.u[end]
        )
        local s2 = solve(NoiseProblem(ou_bridge, (0.0, t_max)); dt = 0.125)
        push_stats!(stats_solver, s2.u)
    end

    @test all(isapprox.(stats_mean(stats_forward), stats_mean(stats_interpolated); atol = 1.0e-3))
    @test all(isapprox.(stats_mean(stats_forward), stats_mean(stats_solver); atol = 1.0e-3))

    @test all(isapprox.(stats_std(stats_forward), stats_std(stats_interpolated); atol = 1.0e-3))
    @test all(isapprox.(stats_std(stats_forward), stats_std(stats_solver); atol = 1.0e-3))
end

@testset "Inplace OU-Bridge" begin
    stats_forward = RunningStats()
    dt = 0.125
    t_max = 20.0
    st = []
    Random.seed!(1234)
    ou = OrnsteinUhlenbeckProcess!([1.5], [0.25], [0.2], 0.0, [0.0])
    for i in 1:500_000
        local prob = NoiseProblem(ou, (0.0, t_max))
        local s = solve(prob; dt = dt)
        push_stats!(stats_forward, s.u)
        if i == 1
            st = s.t
        end
    end

    stats_interpolated = RunningStats()
    Random.seed!(1234)
    ou = OrnsteinUhlenbeckProcess!([1.5], [0.25], [0.2], 0.0, [0.0])
    for i in 1:500_000
        local prob = NoiseProblem(ou, (0.0, t_max))
        local s2 = solve(prob; dt = t_max)
        s2.(st)
        push_stats!(stats_interpolated, s2.u)
    end

    Random.seed!(1234)
    ou = OrnsteinUhlenbeckProcess!([1.5], [0.25], [0.2], 0.0, [0.0])
    prob = NoiseProblem(ou, (0.0, t_max))
    stats_solver = RunningStats()
    for i in 1:500_000
        local s = solve(prob; dt = t_max)
        local ou_bridge = OrnsteinUhlenbeckBridge!(
            [1.5],
            [0.25],
            [0.2],
            0.0,
            t_max,
            [0.0],
            s.u[end]
        )
        local s2 = solve(NoiseProblem(ou_bridge, (0.0, t_max)); dt = 0.125)
        push_stats!(stats_solver, s2.u)
    end

    @test all(isapprox.(stats_mean(stats_forward), stats_mean(stats_interpolated); atol = 1.0e-3))
    @test all(isapprox.(stats_mean(stats_forward), stats_mean(stats_solver); atol = 1.0e-3))

    @test all(isapprox.(stats_std(stats_forward), stats_std(stats_interpolated); atol = 1.0e-3))
    @test all(isapprox.(stats_std(stats_forward), stats_std(stats_solver); atol = 1.0e-3))
end
