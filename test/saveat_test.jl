@testset "saveat Keyword" begin
    using DiffEqNoiseProcess, DiffEqBase, Test, Random

    @testset "Basic saveat functionality with WienerProcess" begin
        # Test that saveat returns only values at specified times
        W = WienerProcess(0.0, 0.0)
        prob = NoiseProblem(W, (0.0, 1.0))
        saveat_times = 0.0:0.1:1.0

        sol = solve(prob; dt = 0.01, saveat = saveat_times)

        # Check that we get exactly the saveat times
        @test length(sol.t) == length(saveat_times)
        @test sol.t ≈ collect(saveat_times)
    end

    @testset "saveat with array input" begin
        W = WienerProcess(0.0, 0.0)
        prob = NoiseProblem(W, (0.0, 1.0))
        saveat_times = [0.0, 0.25, 0.5, 0.75, 1.0]

        sol = solve(prob; dt = 0.01, saveat = saveat_times)

        @test length(sol.t) == 5
        @test sol.t ≈ saveat_times
    end

    @testset "saveat reduces memory compared to save_everystep" begin
        # Solve with fine dt and save_everystep
        W1 = WienerProcess(0.0, 0.0)
        prob1 = NoiseProblem(W1, (0.0, 1.0))
        sol_all = solve(prob1; dt = 0.001)

        # Solve with same dt but saveat at coarser times
        W2 = WienerProcess(0.0, 0.0)
        prob2 = NoiseProblem(W2, (0.0, 1.0))
        sol_saveat = solve(prob2; dt = 0.001, saveat = 0.0:0.1:1.0)

        # saveat solution should have much fewer stored values
        @test length(sol_saveat.t) < length(sol_all.t)
        @test length(sol_saveat.t) == 11  # 0.0, 0.1, ..., 1.0
    end

    @testset "saveat values are consistent with interpolation" begin
        # First solve without saveat
        W = WienerProcess(0.0, 0.0)
        prob = NoiseProblem(W, (0.0, 1.0); seed = 42)
        sol_full = solve(prob; dt = 0.01)

        # Interpolate at saveat times
        saveat_times = [0.0, 0.25, 0.5, 0.75, 1.0]
        interpolated_values = [sol_full(t)[1] for t in saveat_times]

        # Now solve with saveat using same seed
        W2 = WienerProcess(0.0, 0.0)
        prob2 = NoiseProblem(W2, (0.0, 1.0); seed = 42)
        sol_saveat = solve(prob2; dt = 0.01, saveat = saveat_times)

        # Values should match
        for (i, t) in enumerate(saveat_times)
            @test sol_saveat.W[i] ≈ interpolated_values[i]
        end
    end

    @testset "saveat with RealWienerProcess" begin
        W = RealWienerProcess(0.0, 0.0)
        prob = NoiseProblem(W, (0.0, 1.0))
        saveat_times = 0.0:0.2:1.0

        sol = solve(prob; dt = 0.01, saveat = saveat_times)

        @test length(sol.t) == length(saveat_times)
        @test sol.t ≈ collect(saveat_times)
    end

    @testset "saveat with multidimensional WienerProcess" begin
        W = WienerProcess(0.0, [0.0, 0.0])
        prob = NoiseProblem(W, (0.0, 1.0))
        saveat_times = 0.0:0.25:1.0

        sol = solve(prob; dt = 0.01, saveat = saveat_times)

        @test length(sol.t) == 5
        @test all(length(w) == 2 for w in sol.W)
    end

    @testset "saveat with Z process" begin
        W = WienerProcess(0.0, 0.0, 0.0)  # With Z0
        prob = NoiseProblem(W, (0.0, 1.0))
        saveat_times = 0.0:0.25:1.0

        sol = solve(prob; dt = 0.01, saveat = saveat_times)

        @test length(sol.t) == 5
        @test length(sol.Z) == 5
    end

    @testset "saveat skips times outside simulation range" begin
        W = WienerProcess(0.0, 0.0)
        prob = NoiseProblem(W, (0.0, 1.0))
        # Include times outside the simulation range
        saveat_times = [-0.5, 0.0, 0.5, 1.0, 1.5]

        sol = solve(prob; dt = 0.01, saveat = saveat_times)

        # Should only include times within [0.0, 1.0]
        @test length(sol.t) == 3
        @test sol.t ≈ [0.0, 0.5, 1.0]
    end

    @testset "saveat with SimpleWienerProcess" begin
        W = SimpleWienerProcess(0.0, 0.0)
        prob = NoiseProblem(W, (0.0, 1.0))
        saveat_times = 0.0:0.2:1.0

        sol = solve(prob; dt = 0.01, saveat = saveat_times)

        @test length(sol.t) == length(saveat_times)
        @test sol.t ≈ collect(saveat_times)
    end

    @testset "saveat combined with save_everystep=false" begin
        # When save_everystep=false, saveat should still work
        W = WienerProcess(0.0, 0.0; save_everystep = false)
        prob = NoiseProblem(W, (0.0, 1.0))
        saveat_times = 0.0:0.25:1.0

        # This should work even with save_everystep=false because
        # saveat is applied after the simulation
        sol = solve(prob; dt = 0.01, saveat = saveat_times)

        @test length(sol.t) == 5
    end

    @testset "saveat with CorrelatedWienerProcess" begin
        Γ = [1.0 0.5; 0.5 1.0]
        W = CorrelatedWienerProcess(Γ, 0.0, [0.0, 0.0])
        prob = NoiseProblem(W, (0.0, 1.0))
        saveat_times = 0.0:0.25:1.0

        sol = solve(prob; dt = 0.01, saveat = saveat_times)

        @test length(sol.t) == 5
        @test all(length(w) == 2 for w in sol.W)
    end

    @testset "saveat preserves noise process structure" begin
        W = WienerProcess(0.0, 0.0)
        prob = NoiseProblem(W, (0.0, 1.0))
        saveat_times = 0.0:0.25:1.0

        sol = solve(prob; dt = 0.01, saveat = saveat_times)

        # sol should still be a proper noise process
        @test sol isa DiffEqNoiseProcess.NoiseProcess
        # u should be aliased to W
        @test sol.u === sol.W
    end
end
