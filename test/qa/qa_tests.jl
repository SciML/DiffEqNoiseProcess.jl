using Test
using DiffEqNoiseProcess
using JET
using Random

@testset "JET static analysis" begin
    # Note: JET.report_package is skipped here because it creates a virtualized
    # module context that can interfere with subsequent tests in the same process.
    # The optimization tests below provide targeted type stability checks for key functions.

    @testset "Key function optimization" begin
        # Create test instances
        W = WienerProcess(0.0, 0.0)

        # Test key runtime functions for type stability
        @testset "accept_step!" begin
            r = JET.report_opt(
                DiffEqNoiseProcess.accept_step!,
                (typeof(W), Float64, Nothing, Nothing, Bool)
            )
            @test length(JET.get_reports(r)) == 0
        end

        @testset "setup_next_step!" begin
            r = JET.report_opt(
                DiffEqNoiseProcess.setup_next_step!,
                (typeof(W), Nothing, Nothing)
            )
            @test length(JET.get_reports(r)) == 0
        end

        @testset "calculate_step!" begin
            r = JET.report_opt(
                DiffEqNoiseProcess.calculate_step!,
                (typeof(W), Float64, Nothing, Nothing)
            )
            @test length(JET.get_reports(r)) == 0
        end

        @testset "WHITE_NOISE_DIST" begin
            r = JET.report_opt(
                DiffEqNoiseProcess.WHITE_NOISE_DIST,
                (Float64, typeof(W), Float64, Nothing, Nothing, Float64, Random.Xoshiro)
            )
            @test length(JET.get_reports(r)) == 0
        end

        @testset "wiener_randn scalar" begin
            r = JET.report_opt(
                DiffEqNoiseProcess.wiener_randn,
                (Random.Xoshiro, Type{Float64})
            )
            @test length(JET.get_reports(r)) == 0
        end

        @testset "wiener_randn array" begin
            r = JET.report_opt(
                DiffEqNoiseProcess.wiener_randn,
                (Random.Xoshiro, Vector{Float64})
            )
            @test length(JET.get_reports(r)) == 0
        end
    end
end
