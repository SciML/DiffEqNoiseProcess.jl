@testset "Virtual Brownian Tree tests" begin
    using DiffEqNoiseProcess, DiffEqBase, Test, Random
    import RandomNumbers, Random123

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
    rng = RandomNumbers.Random123.Threefry4x() # instantiate PRNG
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
end
