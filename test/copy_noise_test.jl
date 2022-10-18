@testset "Copy Noise test" begin
    using DiffEqNoiseProcess, StochasticDiffEq

    function Base.:(==)(W1::T, W2::T) where {T<:DiffEqNoiseProcess.AbstractNoiseProcess}
        all(getfield(W1, x) == getfield(W2, x) for x in fieldnames(typeof(W1)))
    end

    function Base.:(==)(W1::T, W2::T) where {T<:DiffEqNoiseProcess.NoiseApproximation}
        all(getfield(W1, x) == getfield(W2, x) for x in fieldnames(typeof(W1)) if x ∉ (:source1,:source2))
    end

    for W in (WienerProcess(0.0, 0.0),
              SimpleWienerProcess(0.0, 0.0),
              RealWienerProcess(0.0, 0.0),
              CorrelatedWienerProcess([1.0 0.3; 0.3 1.0],0.0, 0.0),
              GeometricBrownianMotionProcess(0.5, 0.1, 0.0,
                                             1.0),
              OrnsteinUhlenbeckProcess(1.0, 0.2, 1.3, 0.0,
                                       1.0),
              BrownianBridge(0.0, 1.0, 0.0, 1.0))
        W2 = copy(W)
        @test W2 == W
        @test W2 !== W
        @test W2.W === W2.u !== W.W === W.u
    end

    for W in (NoiseFunction(0.0, (u, p, t) -> exp(t)),
              NoiseGrid(0:0.01:1, sin.(0:0.01:1)),
              NoiseWrapper(solve(NoiseProblem(WienerProcess(0.0, 0.0), (0.0, 0.1)),dt=1/10)),
              NoiseApproximation(init(SDEProblem((u, p, t) -> 1.5u, (u, p, t) -> 0.2u, 1.0, (0.0, Inf)), EM(), dt=1/10)),
              VirtualBrownianTree(0.0, 0.0; tree_depth = 3, search_depth = 5),
              BoxWedgeTail(0.0, zeros(2), box_grouping = :Columns)
              )
        W2 = copy(W)
        @test W2 == W
        @test W2 !== W
    end
end
