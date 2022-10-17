@testset "Copy Noises" begin
    using DiffEqNoiseProcess, DiffEqBase

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
              NoiseWrapper(solve(NoiseProblem(WienerProcess(0.0, 0.0), (0.0, 0.1)),dt=1/10))
              )
        W2 = copy(W)
        @test W2 == W
        @test W2 !== W
    end
end
