@testset "NoiseGrid" begin
    using DiffEqNoiseProcess, DiffEqBase, Test

    t = 0:0.001:1
    grid = exp.(t)
    W = NoiseGrid(t, grid)

    dt = 0.1
    calculate_step!(W, dt, nothing, nothing)

    for i in 1:10
        accept_step!(W, dt, nothing, nothing)
    end

    W = NoiseGrid(t, grid)
    prob = NoiseProblem(W, (0.0, 1.0))
    sol = solve(prob; dt = 0.1)

    @test !sol.step_setup
    @test_throws ErrorException accept_step!(sol, dt, nothing, nothing)

    t = 0:0.001:1
    grid = [[exp.(ti) for i in 1:8] for ti in t]
    W = NoiseGrid(t, grid)
    prob = NoiseProblem(W, (0.0, 1.0))
    sol = solve(prob; dt = 0.1)

    dt = 0.001
    t = 0:dt:1
    brownian_values = cumsum([0; [sqrt(dt) * randn() for i in 1:(length(t) - 1)]])
    W = NoiseGrid(t, brownian_values)

    dt = 0.001
    t = 0:dt:1
    brownian_values2 = cumsum(
        [
            [zeros(8)];
            [sqrt(dt) * randn(8) for i in 1:(length(t) - 1)]
        ]
    )
    W = NoiseGrid(t, brownian_values2)
    prob = NoiseProblem(W, (0.0, 1.0))
    sol = solve(prob; dt = 0.1)

    dt = 1 // 1000
    t = 0:dt:1
    W = NoiseGrid(t, brownian_values)
    prob_rational = NoiseProblem(W, (0, 1))
    sol = solve(prob_rational; dt = 1 // 10)
end

@testset "NoiseGrid copies share the grid and not the step state" begin
    using DiffEqNoiseProcess, StochasticDiffEq, Random, Test

    function brownian_grid(n; inplace = false)
        t = collect(range(0.0, 1.0; length = n + 1))
        w = [0.0; cumsum(sqrt(1 / n) .* randn(MersenneTwister(1), n))]
        return inplace ? NoiseGrid(t, [[x, 2x] for x in w]) : NoiseGrid(t, w)
    end

    for inplace in (false, true)
        W = brownian_grid(2^10; inplace)
        W2 = copy(W)
        @test W2.t === W.t
        @test W2.W === W.W
        @test W2.u === W2.W
        @test W2.cur_time !== W.cur_time
        inplace && @test W2.curW !== W.curW && W2.dW !== W.dW

        calculate_step!(W2, 0.1, nothing, nothing)
        accept_step!(W2, 0.1, nothing, nothing)
        @test W.curt == 0.0
        @test W.cur_time[] == 1
        @test W.curW == W.W[1]

        W3 = copy(W)
        copy!(W3, W2)
        @test W3.curt == W2.curt
        @test W3.curW == W2.curW
        @test W3.cur_time[] == W2.cur_time[]
        @test W3.cur_time !== W2.cur_time
        inplace && @test W3.curW !== W2.curW

        small = brownian_grid(2^10; inplace)
        big = brownian_grid(2^16; inplace)
        copy(small)
        copy(big)
        @test (@allocated copy(big)) == (@allocated copy(small))

        f = inplace ? ((du, u, p, t, W) -> (du .= -u .* cos.(5 .* W))) :
            ((u, p, t, W) -> -u * cos(5W))
        u0 = inplace ? [1.0, 1.0] : 1.0
        solve_on(W) = solve(
            RODEProblem{inplace}(f, u0, (0.0, 1.0), noise = W), RandomEM(), dt = 1 / 16
        )
        solve_on(small)
        solve_on(big)
        @test (@allocated solve_on(big)) < 2 * (@allocated solve_on(small))

        prob = RODEProblem{inplace}(f, u0, (0.0, 1.0), noise = brownian_grid(2^10; inplace))
        @test solve(prob, RandomEM(), dt = 1 / 16).u == solve(prob, RandomEM(), dt = 1 / 16).u
    end
end
