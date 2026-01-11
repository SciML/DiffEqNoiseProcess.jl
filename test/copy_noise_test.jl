@testset "Copy Noise test" begin
    using DiffEqNoiseProcess, StochasticDiffEq
    using Optim  # Required for BoxWedgeTail with sqeezing=true

    # define a temporary equality suitable for comparing noise types (tab-completed by `\boxminus<tab>`)
    ⊟(W1, W2) = (W1 == W2)
    function ⊟(W1::T, W2::T) where {T <: DiffEqNoiseProcess.AbstractNoiseProcess}
        equal = true
        for x in fieldnames(T)
            xequal = true
            if getfield(W2, x) isa DiffEqNoiseProcess.ResettableStacks.ResettableStack
                xequal &= all(
                    getfield(getfield(W1, x), y) ⊟ getfield(getfield(W2, x), y)
                        for y in (:cur, :numResets, :data)
                )
            elseif getfield(W2, x) isa DiffEqNoiseProcess.RSWM
                xequal &= all(
                    getfield(getfield(W1, x), y) ⊟ getfield(getfield(W2, x), y)
                        for y in (:discard_length, :adaptivealg)
                )
            elseif !ismutable(getfield(W1, x)) || getfield(W1, x) isa AbstractArray ||
                    getfield(W1, x) === nothing || getfield(W2, x) === nothing
                xequal &= (getfield(W1, x) ⊟ getfield(W2, x))
            end
            if xequal != true
                @info "$x::$(typeof(getfield(W1, x))) with value W1.$x = $(getfield(W1, x)) and W2.$x = $(getfield(W2, x)) in $(first(split(string(T), '}')))"
            end
            equal &= xequal
        end
        equal
    end

    i = 0
    for W in (
            WienerProcess(0.0, 0.0),
            SimpleWienerProcess(0.0, 0.0),
            RealWienerProcess(0.0, 0.0),
            CorrelatedWienerProcess([1.0 0.3; 0.3 1.0], 0.0, 0.0),
            GeometricBrownianMotionProcess(
                0.5, 0.1, 0.0,
                1.0
            ),
            OrnsteinUhlenbeckProcess(
                1.0, 0.2, 1.3, 0.0,
                1.0
            ),
            BrownianBridge(0.0, 1.0, 0.0, 1.0),
        )
        W2 = deepcopy(W)
        @test typeof(W2) == typeof(W)
        copy!(W2, W)
        @test W2 ⊟ W
        @test copy(W) ⊟ W
        @test W2 !== W
        @test W2.W === W2.u !== W.W === W.u
    end

    for (W1, W2) in (
            (WienerProcess(0.0, 0.0), WienerProcess(1.0, 1.0)),
            (SimpleWienerProcess(0.0, 0.0), SimpleWienerProcess(1.0, 1.0)),
            (RealWienerProcess(0.0, 0.0), RealWienerProcess(1.0, 1.0)),
        )
        W = deepcopy(W1)
        @test typeof(W2) == typeof(W1)
        @test W ⊟ W1
        copy!(W2, W1)
        @test W2 ⊟ W1
        @test copy(W1) ⊟ W1
        @test W2 !== W1
        @test W2.W === W2.u !== W1.W === W1.u
    end

    for W in (
            NoiseFunction(0.0, (u, p, t) -> exp(t)),
            NoiseTransport(0.0, (u, p, t, Y) -> exp(t), (rng) -> nothing),
            NoiseGrid(0:0.01:1, sin.(0:0.01:1)),
            NoiseWrapper(
                solve(
                    NoiseProblem(WienerProcess(0.0, 0.0), (0.0, 0.1)),
                    dt = 1 / 10
                )
            ),
            NoiseApproximation(
                init(
                    SDEProblem(
                        (u, p, t) -> 1.5u, (u, p, t) -> 0.2u, 1.0,
                        (0.0, Inf)
                    ),
                    EM(), dt = 1 / 10
                )
            ),
            VirtualBrownianTree(0.0, 0.0; tree_depth = 3, search_depth = 5),
            BoxWedgeTail(0.0, zeros(2), box_grouping = :Columns),
        )
        W2 = deepcopy(W)
        @test typeof(W2) == typeof(W)
        copy!(W2, W)
        @test W2 ⊟ W
        @test copy(W) ⊟ W
        @test W2 !== W
    end
end
