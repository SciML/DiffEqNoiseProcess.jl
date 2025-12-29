# Regression test for issue #231 / PR #230
# Tests that SimpleWienerProcess correctly passes dW as first argument to dist/bridge functions
# Before PR #230, this would fail with BigFloat because WHITE_NOISE_DIST uses dW to determine output type

@testset "SimpleNoiseProcess BigFloat" begin
    using DiffEqNoiseProcess

    # Test 1: Out-of-place SimpleWienerProcess with BigFloat scalar
    @testset "Scalar BigFloat (out-of-place)" begin
        W = SimpleWienerProcess(big(0.0), big(0.0))
        dt = big(0.1)

        # This should work - calculate_step! needs to pass dW as first argument
        # Regression test for PR #230
        calculate_step!(W, dt, nothing, nothing)
        @test W.dW isa BigFloat

        # Accept step should also work
        accept_step!(W, dt, nothing, nothing)
        @test W.curW isa BigFloat
    end

    # Test 2: Out-of-place SimpleWienerProcess with BigFloat array
    @testset "Array BigFloat (out-of-place)" begin
        W = SimpleWienerProcess(big(0.0), zeros(BigFloat, 3))
        dt = big(0.1)

        calculate_step!(W, dt, nothing, nothing)
        @test eltype(W.dW) == BigFloat

        accept_step!(W, dt, nothing, nothing)
        @test eltype(W.curW) == BigFloat
    end

    # Test 3: In-place SimpleWienerProcess with BigFloat array
    @testset "Array BigFloat (in-place)" begin
        W = SimpleWienerProcess!(big(0.0), zeros(BigFloat, 3))
        dt = big(0.1)

        calculate_step!(W, dt, nothing, nothing)
        @test eltype(W.dW) == BigFloat

        accept_step!(W, dt, nothing, nothing)
        @test eltype(W.curW) == BigFloat
    end
end
