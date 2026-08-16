@testset "AbstractNoiseProcess lifecycle interface" begin
    using DiffEqNoiseProcess
    using SciMLBase: AbstractNoiseProcess

    source = WienerProcess(0.0, 0.0; save_everystep = true)
    calculate_step!(source, 0.1, nothing, nothing)
    accept_step!(source, 0.1, nothing, nothing, false)

    rejectable = AbstractNoiseProcess[
        WienerProcess(0.0, 0.0),
        NoiseFunction(0.0, (u, p, t) -> t),
        NoiseTransport(0.0, (u, p, t, rv) -> rv * t, _ -> 1.0),
        NoiseGrid(collect(0.0:0.1:1.0), collect(0.0:0.1:1.0)),
        NoiseWrapper(source),
    ]

    for W in rejectable
        @testset "$(nameof(typeof(W)))" begin
            @test W(0.0) isa Tuple
            @test setup_next_step!(W, nothing, nothing) === nothing
            calculate_step!(W, 0.1, nothing, nothing)
            @test reject_step!(W, 0.05, nothing, nothing) === nothing
            @test accept_step!(W, 0.05, nothing, nothing, false) === nothing
            @test save_noise!(W) === nothing
        end
    end

    @testset "SimpleNoiseProcess rejects adaptive stepping" begin
        W = SimpleWienerProcess(0.0, 0.0)
        @test W isa AbstractNoiseProcess
        @test setup_next_step!(W, nothing, nothing) === nothing
        calculate_step!(W, 0.1, nothing, nothing)
        @test_throws ErrorException reject_step!(W, 0.05, nothing, nothing)
        @test accept_step!(W, 0.1, nothing, nothing, false) === nothing
        @test save_noise!(W) === nothing
    end

    @testset "In-place callable contract" begin
        processes = AbstractNoiseProcess[
            WienerProcess!(0.0, zeros(2)),
            NoiseFunction(
                0.0,
                (out, u, p, t) -> fill!(out, t),
                noise_prototype = zeros(2),
            ),
        ]
        for W in processes
            out = zeros(2)
            W(out, nothing, nothing, nothing, 0.0)
            @test out == zeros(2)
        end
    end
end
