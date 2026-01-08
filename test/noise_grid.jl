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

    # Test for issue #136: stepping past domain with dt that doesn't evenly divide the time span
    # When dt doesn't evenly divide (tspan[2] - tspan[1]), floating point errors can cause
    # the solver to think it needs to step past the domain boundary
    @testset "Issue #136 - stepping past domain" begin
        # Original issue case: dt=0.01 with grid step 0.1
        tgrid = 0.0:0.1:10.0
        brownian_noise = [[0.0, 0.0] for _ in 1:length(tgrid)]
        W = NoiseGrid(collect(tgrid), brownian_noise)
        prob = NoiseProblem(W, (tgrid[begin], tgrid[end]))
        sol = solve(prob; dt = 0.01)
        @test sol.curt == 10.0

        # Case where dt does not evenly divide the time span
        tgrid2 = 0.0:0.1:1.0
        brownian_noise2 = [[0.0, 0.0] for _ in 1:length(tgrid2)]
        W2 = NoiseGrid(collect(tgrid2), brownian_noise2)
        prob2 = NoiseProblem(W2, (0.0, 1.0))
        sol2 = solve(prob2; dt = 0.03)
        @test sol2.curt == 1.0

        # Case with scalar noise
        tgrid3 = 0.0:0.1:10.0
        brownian_noise3 = [0.0 for _ in 1:length(tgrid3)]
        W3 = NoiseGrid(collect(tgrid3), brownian_noise3)
        prob3 = NoiseProblem(W3, (tgrid3[begin], tgrid3[end]))
        sol3 = solve(prob3; dt = 0.01)
        @test sol3.curt == 10.0
    end
end
