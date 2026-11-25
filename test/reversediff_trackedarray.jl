using DiffEqNoiseProcess
using StochasticDiffEq
using SciMLBase
using SciMLSensitivity
using ReverseDiff
using ForwardDiff
using Zygote
using Random, Test

@testset "ReverseDiff tracked noise storage" begin
    # A TrackedArray prototype cannot back an in-place noise process: ReverseDiff
    # >= 1.18.6 forbids all mutation of TrackedArray, while setup_next_step! and
    # friends fill!/copyto!/broadcast into the W0-derived buffers.
    proto = ReverseDiff.track(zeros(2))
    W = WienerProcess!(0.0, proto)
    @test W.dW isa Array{<:ReverseDiff.TrackedReal}
    @test W.curW isa Array{<:ReverseDiff.TrackedReal}
    @test fill!(W.dW, zero(eltype(W.dW))) === W.dW

    function seeded_wiener_step(w0, z0)
        W = WienerProcess!(0.0, w0, z0; rng = Xoshiro(12345))
        calculate_step!(W, 0.25, nothing, nothing)
        return W.dW, W.dZ
    end
    plain_dW, plain_dZ = seeded_wiener_step(zeros(16), zeros(16))
    tracked_dW, tracked_dZ = seeded_wiener_step(
        ReverseDiff.track(zeros(16)), ReverseDiff.track(zeros(16))
    )
    @test map(x -> getfield(x, :value), tracked_dW) == plain_dW
    @test map(x -> getfield(x, :value), tracked_dZ) == plain_dZ

    aux_tracked_dW, aux_tracked_dZ = seeded_wiener_step(
        zeros(16), ReverseDiff.track(zeros(16))
    )
    @test aux_tracked_dW == plain_dW
    @test map(x -> getfield(x, :value), aux_tracked_dZ) == plain_dZ

    # Out-of-place processes keep TrackedArray storage: they rebind, not mutate,
    # and wiener_randn on a TrackedArray prototype draws through the tape.
    Woop = WienerProcess(0.0, proto)
    @test Woop.dW isa ReverseDiff.TrackedArray

    auxiliary_sum(w0, z0) = begin
        W = WienerProcess!(0.0, w0, z0; rng = MersenneTwister(1))
        calculate_step!(W, 0.1, nothing, nothing)
        accept_step!(W, 0.1, nothing, nothing)
        sum(W.Z[end])
    end
    @test ReverseDiff.gradient(z0 -> auxiliary_sum([1.0, -2.0], z0), [3.0, -4.0]) ==
        [1.0, 1.0]

    both_tracked = WienerProcess!(
        0.0, ReverseDiff.track([1.0, -2.0]), ReverseDiff.track([3.0, -4.0]);
        rng = MersenneTwister(1)
    )
    @test both_tracked.dW isa Array{<:ReverseDiff.TrackedReal}
    @test both_tracked.dZ isa Array{<:ReverseDiff.TrackedReal}
    calculate_step!(both_tracked, 0.1, nothing, nothing)
    accept_step!(both_tracked, 0.1, nothing, nothing)

    function rode_f(du, u, p, t, W)
        du[1] = p[1] * u[1] * sin(W[1] - W[2])
        du[2] = p[2] * u[2] * cos(W[1] + W[2])
        return nothing
    end
    seed = 12345
    u0 = [1.0, 1.0]
    p = [2.0, -2.0]
    tspan = (0.0, 0.5)
    dt = 1.0e-3
    t = tspan[1]:0.1:tspan[2]
    prob = RODEProblem(rode_f, u0, tspan, p)
    Random.seed!(seed)
    du0, dp = Zygote.gradient(
        (u0, p) -> sum(
            Array(solve(prob, RandomEM(); dt, u0, p, saveat = t, sensealg = ReverseDiffAdjoint())) .^ 2 / 2
        ),
        u0, p
    )
    Random.seed!(seed)
    df = ForwardDiff.gradient(
        θ -> sum(Array(solve(prob, RandomEM(); dt, u0 = θ[1:2], p = θ[3:4], saveat = t)) .^ 2 / 2),
        [u0; p]
    )
    @info du0, dp, df
    @test du0 ≈ df[1:2] rtol = 1.0e-7
    @test dp ≈ df[3:4] rtol = 1.0e-7

    function wide_rode_f(du, u, p, t, W)
        @. du = p[1] * u * sin(W) + p[2] * cos(W)
        return nothing
    end
    wide_u0 = collect(range(0.8, 1.1; length = 16))
    wide_p = [1.2, -0.8]
    wide_t = 0.0:0.1:0.5
    wide_prob = RODEProblem(wide_rode_f, wide_u0, (0.0, 0.5), wide_p)
    wide_loss(p; sensealg = nothing) = sum(
        abs2, Array(
            solve(
                wide_prob, RandomEM(); dt = 0.01, p, saveat = wide_t,
                seed = UInt64(12345), sensealg
            )
        )
    ) / 2
    wide_adjoint = Zygote.gradient(
        p -> wide_loss(p; sensealg = ReverseDiffAdjoint()), wide_p
    )[1]
    fd_gradient(h) = [
        (wide_loss(wide_p + h .* [1.0, 0.0]) - wide_loss(wide_p - h .* [1.0, 0.0])) / (2h),
        (wide_loss(wide_p + h .* [0.0, 1.0]) - wide_loss(wide_p - h .* [0.0, 1.0])) / (2h),
    ]
    fd_h = fd_gradient(1.0e-5)
    fd_h2 = fd_gradient(0.5e-5)
    @test fd_h ≈ fd_h2 rtol = 1.0e-5
    @test wide_adjoint ≈ fd_h2 rtol = 1.0e-4
end
