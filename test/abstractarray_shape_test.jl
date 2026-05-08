@testset "AbstractNoiseProcess T/N parameterization" begin
    using DiffEqNoiseProcess, StochasticDiffEq, Test

    # AbstractNoiseProcess <: AbstractDiffEqArray <: AbstractVectorOfArray <: AbstractArray{T, N}
    # so eltype must be the inner scalar (Float64) and ndims must be
    # `1 + ndims(inner_sample)` (one inner dim per state component, plus time).

    @testset "scalar noise: T=Float64, N=1" begin
        ts = collect(0.0:0.1:1.0)
        Ws = randn(length(ts))

        @testset "NoiseGrid" begin
            W = NoiseGrid(ts, Ws)
            @test eltype(W) === Float64
            @test ndims(W) == 1
        end

        @testset "NoiseProcess" begin
            W = WienerProcess(0.0, 0.0, 0.0)
            @test eltype(W) === Float64
            @test ndims(W) == 1
        end

        @testset "SimpleNoiseProcess" begin
            W = SimpleWienerProcess(0.0, 0.0, 0.0)
            @test eltype(W) === Float64
            @test ndims(W) == 1
        end

        @testset "NoiseFunction" begin
            W = NoiseFunction(0.0, (u, p, t) -> exp(t))
            @test eltype(W) === Float64
            @test ndims(W) == 1
        end

        @testset "NoiseTransport" begin
            W = NoiseTransport(0.0, (u, p, t, Y) -> exp(t), (rng) -> nothing)
            @test eltype(W) === Float64
            @test ndims(W) == 1
        end

        @testset "VirtualBrownianTree" begin
            W = VirtualBrownianTree(0.0, 0.0; tree_depth = 3, search_depth = 5)
            @test eltype(W) === Float64
            @test ndims(W) == 1
        end

        @testset "NoiseApproximation" begin
            integ = init(
                SDEProblem(
                    (u, p, t) -> 1.5u, (u, p, t) -> 0.2u, 1.0,
                    (0.0, Inf)
                ),
                EM(); dt = 1 / 10
            )
            W = NoiseApproximation(integ)
            @test eltype(W) === Float64
            @test ndims(W) == 1
        end
    end

    @testset "vector noise: T=Float64, N=2" begin
        ts = collect(0.0:0.1:1.0)
        Ws = [randn(2) for _ in ts]

        @testset "NoiseGrid" begin
            W = NoiseGrid(ts, Ws)
            @test eltype(W) === Float64
            @test ndims(W) == 2
            @test size(W) == (2, length(ts))
            # Indexed access returns the underlying scalar.
            @test W[1, 1] === Ws[1][1]
            @test W[2, end] === Ws[end][2]
        end

        @testset "NoiseProcess" begin
            W = WienerProcess!(0.0, zeros(2), zeros(2))
            @test eltype(W) === Float64
            @test ndims(W) == 2
        end

        @testset "SimpleNoiseProcess" begin
            W = SimpleWienerProcess!(0.0, zeros(2), zeros(2))
            @test eltype(W) === Float64
            @test ndims(W) == 2
        end

        @testset "NoiseFunction" begin
            W = NoiseFunction(0.0, (u, p, t) -> [exp(t), exp(2t)])
            @test eltype(W) === Float64
            @test ndims(W) == 2
        end

        @testset "NoiseTransport" begin
            W = NoiseTransport(
                0.0, (u, p, t, Y) -> Y .* exp(t), (rng) -> randn(rng, 2)
            )
            @test eltype(W) === Float64
            @test ndims(W) == 2
        end

        @testset "VirtualBrownianTree" begin
            W = VirtualBrownianTree(
                0.0, zeros(2); tree_depth = 3, search_depth = 5
            )
            @test eltype(W) === Float64
            @test ndims(W) == 2
        end

        @testset "NoiseApproximation" begin
            integ = init(
                SDEProblem(
                    (du, u, p, t) -> (du .= 1.5 .* u),
                    (du, u, p, t) -> (du .= 0.2 .* u),
                    [1.0, 1.0], (0.0, Inf)
                ),
                EM(); dt = 1 / 10
            )
            W = NoiseApproximation(integ)
            @test eltype(W) === Float64
            @test ndims(W) == 2
        end
    end

    @testset "matrix noise: T=Float64, N=3" begin
        ts = collect(0.0:0.1:1.0)
        Ws = [randn(2, 3) for _ in ts]

        @testset "NoiseGrid" begin
            W = NoiseGrid(ts, Ws)
            @test eltype(W) === Float64
            @test ndims(W) == 3
            @test size(W) == (2, 3, length(ts))
            @test W[1, 1, 1] === Ws[1][1, 1]
            @test W[2, 3, end] === Ws[end][2, 3]
        end

        @testset "NoiseFunction" begin
            W = NoiseFunction(0.0, (u, p, t) -> ones(2, 3) .* t)
            @test eltype(W) === Float64
            @test ndims(W) == 3
        end
    end
end
