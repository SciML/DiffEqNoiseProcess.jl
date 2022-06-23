@testset "NoiseFunction" begin

    using DiffEqNoiseProcess, DiffEqBase, Test
    f = (u, p, t) -> exp(t)

    W = NoiseFunction(0.0, f)

    dt = 0.1
    calculate_step!(W, dt, nothing, nothing)

    for i = 1:10
        accept_step!(W, dt, nothing, nothing)
    end

    @test W(W.curt + W.dt)[1] - W.curW == W.dW

    prob = NoiseProblem(W, (0.0, 1.0))
    sol = solve(prob; dt = 0.1)

    @test sol(W.curt + W.dt)[1] - W.curW == W.dW

    f = (out, u, p, t) -> (out .= exp(t))
    W = NoiseFunction(0.0, f, noise_prototype = rand(4))
    prob = NoiseProblem(W, (0.0, 1.0))
    sol = solve(prob; dt = 0.1)

    out1 = rand(4)
    out2 = nothing
    sol(out1, out2, nothing, nothing, 0.5)

    @test out1 == sol(0.5)[1]
end
