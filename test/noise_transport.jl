@testset "NoiseTransport" begin
    using DiffEqNoiseProcess, DiffEqBase, Random, Test

    f = (u, p, t, Y) -> exp(Y * t)

    W = NoiseTransport(0.0, f, randn)

    dt = 0.1
    calculate_step!(W, dt, nothing, nothing)

    for i in 1:10
        accept_step!(W, dt, nothing, nothing)
    end

    @test W(W.curt + W.dt)[1] - W.curW == W.dW

    prob = NoiseProblem(W, (0.0, 1.0))
    sol = solve(prob; dt = 0.1)

    @test sol(nothing, nothing, sol.curt + sol.dt, sol.rv)[1] - sol.curW == sol.dW

    f = (out, u, p, t, Y) -> (out .= exp(Y * t))
    W = NoiseTransport(0.0, f, randn, 1.0, noise_prototype = rand(4))
    prob = NoiseProblem(W, (0.0, 1.0))
    sol = solve(prob; dt = 0.1)

    out1 = rand(4)
    out2 = nothing
    sol(out1, out2, nothing, nothing, 0.5, sol.rv)

    @test out1 == sol(0.5)[1] == W(nothing, nothing, 0.5, sol.rv)[1]

end
