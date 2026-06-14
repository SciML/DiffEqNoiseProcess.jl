@safetestset "save_everystep Keyword" begin
    #Test whether the result of the process is dependent on 'save_everystep'.
    using DiffEqNoiseProcess, DiffEqBase, Test, Statistics, Random
    processes = [
        OrnsteinUhlenbeckProcess(1.0, 1.0, 0.3, 0.0, 0.0, nothing),
        WienerProcess(0.0, 0.0, nothing),
        CorrelatedWienerProcess([1.0 0.0; 0.0 1.0], 0.0, [0.0; 0.0], nothing),
        GeometricBrownianMotionProcess(1.0, 1.0, 0.0, 0.0, nothing),
    ]

    @testset "Noise_process = $(proc.dist)" for proc in processes
        cproc = deepcopy(proc)
        cproc.save_everystep = true
        prob = NoiseProblem(cproc, (0.0, 1.0); seed = 1234)
        sol_save = solve(prob; dt = 0.1)

        cproc = deepcopy(proc)
        cproc.save_everystep = false
        prob = NoiseProblem(cproc, (0.0, 1.0); seed = 1234)
        sol_nosave = solve(prob; dt = 0.1)

        @test sol_save.curW == sol_nosave.curW
    end
end

@safetestset "Noise solution lands exactly on the endpoint" begin
    # A NoiseProblem solve must stop exactly on tspan[2] and never step past it.
    # Floating point drift accumulated over many steps is ~eps(tspan[2]), which can be
    # far larger than eps(dt) for small dt, so the end correction must be scaled by the
    # magnitude of the time values. Otherwise the solve overshoots the requested endpoint.
    using DiffEqNoiseProcess, DiffEqBase, Test, Random
    processes = [
        OrnsteinUhlenbeckProcess(1.0, 1.0, 0.3, 0.0, 0.0, nothing),
        WienerProcess(0.0, 0.0, nothing),
        CorrelatedWienerProcess([1.0 0.0; 0.0 1.0], 0.0, [0.0; 0.0], nothing),
        GeometricBrownianMotionProcess(1.0, 1.0, 0.0, 0.0, nothing),
    ]
    # tspan/dt combinations chosen so that dt does not divide the span cleanly and the
    # accumulated drift is large compared to eps(dt) (the previously buggy regime).
    cases = [
        ((0.0, 1.0), 1.0e-4),
        ((0.0, 1.0), 0.1),
        ((0.0, 2.3), 0.01),
        ((0.0, 1.0), 1 / 3),
    ]
    @testset "process = $(proc.dist)" for proc in processes
        @testset "tspan = $tspan, dt = $dt" for (tspan, dt) in cases
            cproc = deepcopy(proc)
            prob = NoiseProblem(cproc, tspan; seed = 1234)
            sol = solve(prob; dt = dt)
            @test sol.t[end] == tspan[2]
            @test sol.curt == tspan[2]
            @test sol.t[end] <= tspan[2]
            @test issorted(sol.t)
        end
    end
end
