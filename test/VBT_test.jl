@testset "Virtual Brownian Tree tests" begin
    using DiffEqNoiseProcess, DiffEqBase, StochasticDiffEq, Test, Random
    import Random123

    W = VirtualBrownianTree(0.0, 0.0; tree_depth = 3, search_depth = 5)
    @test isinplace(W) == false
    @test length(W.seeds) == 2^(W.tree_depth)
    @test length(W.W) == 2^(W.tree_depth) + 1
    @test length(W.t) == 2^(W.tree_depth) + 1
    @test length(unique(W.seeds)) == length(W.seeds)

    # step to cached value
    dt = W.t[2] - W.t[1]
    calculate_step!(W, dt, nothing, nothing)
    @test W.dW == W.W[2] - W.W[1]

    for i in 1:(length(W.t) - 1)
        accept_step!(W, dt, nothing, nothing)
    end
    @test_throws ErrorException accept_step!(W, dt, nothing, nothing)

    # test interpolation/binary search
    rng = Random123.Threefry4x() # instantiate PRNG
    rngcopy = copy(rng)
    W = VirtualBrownianTree(0.0, 0.0; tree_depth = 1, search_depth = 5, rng = rng)
    t = (W.t[2] - W.t[1]) / 8
    # VBT but with deeper tree such that no binary search is required
    Wcopy = VirtualBrownianTree(0.0, 0.0; tree_depth = 4, search_depth = 5, rng = rngcopy)
    @test W(t) == Wcopy(t)

    W = VirtualBrownianTree(0.0, 0.0; tree_depth = 1, rng = rng)
    for i in 1:10
        accept_step!(W, 0.1, nothing, nothing)
    end

    # test NoiseProblem
    W = VirtualBrownianTree(0.0, 0.0; tree_depth = 1, rng = rng)
    prob = NoiseProblem(W, (0.0, 1.0))
    sol = solve(prob; dt = 0.05)

    # test rational numbers
    W = VirtualBrownianTree(0 // 1, 0.0; tree_depth = 1, rng = rng)
    prob = NoiseProblem(W, (0, 1))
    sol = solve(prob; dt = 1 // 2)

    W = VirtualBrownianTree!(0.0, zeros(8); tree_depth = 3, search_depth = 5)
    @test isinplace(W) == true
    dt = W.t[2] - W.t[1]
    calculate_step!(W, dt, nothing, nothing)
    @test W.dW == W.W[2] - W.W[1]

    for i in 1:(length(W.t) - 1)
        accept_step!(W, dt, nothing, nothing)
    end

    # test memory size
    w0 = zeros(1_000_000)
    W = VirtualBrownianTree(0.0, w0; tree_depth = 0)
    # size = 4*W0 to store W0, Wend, curW, and dW
    @test isapprox(Base.summarysize(W) / Base.summarysize(w0), 4, rtol = 1e-4)

    prob = NoiseProblem(W, (0.0, 1.0))
    sol = solve(prob; dt = 1 / 10)

    @test isapprox(Base.summarysize(sol) / Base.summarysize(w0), 4, rtol = 1e-4)

    # test reproducibility
    f = (u, p, t) -> 1.01u
    g = (u, p, t) -> 1.01u
    dt = 1 // 2^4

    ## forward reproducibility
    W = VirtualBrownianTree(0.0, 0.0; tree_depth = 5, tend = 2.0)
    prob = SDEProblem(f, g, 1.0, (0.0, 1.0), noise = W)

    sol1 = solve(prob, EM(), dt = dt, save_noise = true)
    sol2 = solve(prob, EM(), dt = dt)

    @test sol1.W.W ≈ sol2.W.W
    @test sol1.u ≈ sol2.u

    sol2 = solve(prob, EM(), dt = dt)
    @test sol1.W.W ≈ sol2.W.W
    @test sol1.u ≈ sol2.u

    prob2 = SDEProblem(f, g, sol1(0.5), (0.5, 1.0), noise = W)

    sol2 = solve(prob2, EM(), dt = dt, save_noise = true)

    @test sol1.W.W ≈ sol2.W.W
    @test sol2.u[1] == sol1(0.5)
    @test sol2.u[2] ≈ sol1(0.5 + dt)
    @test sol1(sol2.t).u ≈ sol2.u

    ## backward reproducibility
    W = VirtualBrownianTree(-2.0, 0.0; tree_depth = 6, tend = 2.0)
    prob = SDEProblem(f, g, 1.0, (1.0, 0.0), noise = W)

    sol1 = solve(prob, EM(), dt = dt, save_noise = true)
    sol2 = solve(prob, EM(), dt = dt)

    @test sol1.W.W ≈ sol2.W.W
    @test sol1.u ≈ sol2.u

    sol2 = solve(prob, EM(), dt = dt)
    @test sol1.W.W ≈ sol2.W.W
    @test sol1.u ≈ sol2.u

    prob2 = SDEProblem(f, g, sol1(0.5), (0.5, 0.0), noise = W)

    sol2 = solve(prob2, EM(), dt = dt, save_noise = true)

    @test sol1.W.W ≈ sol2.W.W
    @test sol2.u[1] == sol1(0.5)
    @test sol2.u[2] ≈ sol1(0.5 - dt)
    @test sol1(sol2.t).u ≈ sol2.u

    ## forward vs backward reproducibility with reduced diffusion
    f = (u, p, t) -> 1.01u
    g = (u, p, t) -> 0.11u
    W = VirtualBrownianTree(-2.0, 0.0; tree_depth = 6, tend = 2.0)
    dt = 1 / 2^8
    prob1 = SDEProblem(f, g, 1.0, (0.0, 1.0), noise = W)
    sol1 = solve(prob1, EM(), dt = dt, save_noise = true)
    prob2 = SDEProblem(f, g, sol1(1.0), (1.0, 0.0), noise = W)
    sol2 = solve(prob2, EM(), dt = dt, save_noise = true)
    sol3 = solve(prob2, EM(), dt = dt, save_noise = true)
    @test sol1.W.W ≈ sol2.W.W
    @test sol2.u ≈ sol3.u
    @test sol2.u≈reverse(sol1.u) rtol=0.01
end
