@testset "NoiseTransport" begin
    using DiffEqNoiseProcess, Distributions, Random, Test

    # test steps
    f = (u, p, t, Y) -> exp(Y * t)

    W = NoiseTransport(0.0, f, randn)

    dt = 0.1
    calculate_step!(W, dt, nothing, nothing)

    for i in 1:10
        accept_step!(W, dt, nothing, nothing)
    end

    @test W(W.curt + W.dt)[1] - W.curW == W.dW

    # test solving it
    prob = NoiseProblem(W, (0.0, 1.0))
    sol = solve(prob; dt = 0.1)

    @test sol(nothing, nothing, sol.curt + sol.dt, sol.rv)[1] - sol.curW == sol.dW

    # compare with NoiseFunction
    f1 = (u, p, t, Y) -> exp(t)
    f2 = (u, p, t) -> exp(t)
    RV = (rng) -> nothing
    W = NoiseTransport(0.0, f1, RV)
    Y = NoiseFunction(0.0, f2)

    tt = 0.0:0.1:1.0
    @test W.(tt) == Y.(tt)

    # system with scaler noise
    f = (out, u, p, t, Y) -> (out .= exp(Y * t))
    W = NoiseTransport(0.0, f, randn, 1.0, noise_prototype = rand(4))
    prob = NoiseProblem(W, (0.0, 1.0))
    sol = solve(prob; dt = 0.1)

    out1 = rand(4)
    out2 = nothing
    sol(out1, out2, nothing, nothing, 0.5, sol.rv)

    @test out1 == sol(0.5)[1] == W(nothing, nothing, 0.5, sol.rv)[1]

    # scalar with random vector
    f = (u, p, t, rv) -> sin(p[1] * t + rv[1]) + cos(p[2] * t + rv[2])
    t0 = 0.0
    rv = randn(2)
    p = (π, 2π)
    @test_nowarn NoiseTransport(t0, f, randn!, rv, noise_prototype = f(nothing, p, t0, rv))

    # system with random variable from Distributions.jl
    f! = (out, u, p, t, rv) -> (out .= sin.(rv * t))
    RV(rng) = rand(rng, Beta(2, 3))
    @test_nowarn NoiseTransport(t0, f!, RV; noise_prototype = zeros(4))

    # system with random vector
    f! = function (out, u, p, t, v)
        out[1] = sin(v[1] * t)
        out[2] = sin(t + v[2])
        out[3] = cos(t) * v[1] + sin(t) * v[2]
        nothing
    end

    RV! = (rng, v) -> (v[1] = randn(rng); v[2] = rand(rng))
    rv = zeros(2)

    @test_nowarn NoiseTransport(t0, f!, RV!, rv, noise_prototype = zeros(3))
end
