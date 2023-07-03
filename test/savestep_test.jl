@testset "save_everystep Keyword" begin
    #Test whether the result of the process is dependent on 'save_everystep'.
    using DiffEqNoiseProcess, DiffEqBase, Test, Statistics, Random
    processes = [OrnsteinUhlenbeckProcess(1.0, 1.0, 0.3, 0.0, 0.0, nothing),
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
