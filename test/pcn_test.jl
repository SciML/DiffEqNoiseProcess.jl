@testset "preconditioned Crank Nicolson tests" begin
    using DiffEqNoiseProcess, Test, Random
    using Statistics
    using DiffEqBase
    using DiffEqBase.EnsembleAnalysis

    ##
    # Tests of pCN
    ##

    W = WienerProcess(0.0, 0.0, 0.0)

    dt = 0.1
    calculate_step!(W, dt, nothing, nothing)

    for i in 1:10
        accept_step!(W, dt, nothing, nothing)
    end

    _W = deepcopy(W)

    # test with ρ=1
    W2 = pCN!(_W, 1.0)
    W2a = pCN(W, 1.0)
    WWrapper = NoiseWrapper(_W)

    # test source
    @test W2.source.W == WWrapper.source.W
    @test W2.source.Z == WWrapper.source.Z
    @test W2.source.W == W.W
    @test W2.source.u == W.u
    @test W2.source.Z == W.Z
    @test W2.source.W == W2a.source.W
    @test W2.source.Z == W2a.source.Z

    # multi dimensional
    W = WienerProcess(0.0, zeros(4), zeros(4))

    dt = 0.1
    calculate_step!(W, dt, nothing, nothing)

    for i in 1:10
        accept_step!(W, dt, nothing, nothing)
    end

    _W = deepcopy(W)

    W2 = pCN!(_W, 1.0)
    W2a = pCN(W, 1.0)
    WWrapper = NoiseWrapper(_W)

    # test source
    @test W2.source.W == WWrapper.source.W
    @test W2.source.Z == WWrapper.source.Z
    @test W2.source.W == _W.W
    @test W2.source.u == _W.u
    @test W2.source.Z == _W.Z
    @test W2.source.W == W.W
    @test W2.source.Z == W.Z
    @test W2.source.W == W2a.source.W
    @test W2.source.Z == W2a.source.Z

    # test ρ!=0 and ρ!=1
    _W = deepcopy(W)
    W3 = pCN!(_W, 0.2)
    W3a = pCN(W, 0.2)

    # test source

    @test W3.source.W != W.W
    @test W3.source.u != W.u
    @test W3.source.Z == W.Z # no action on auxilary process
    @test W3.source.W == W3a.source.W
    @test W3.source.u == W3a.source.u
    @test W3.source.Z == W3a.source.Z

    # inplace
    W = WienerProcess!(0.0, zeros(4), zeros(4))

    dt = 0.1
    calculate_step!(W, dt, nothing, nothing)

    for i in 1:10
        accept_step!(W, dt, nothing, nothing)
    end

    _W = deepcopy(W)

    # test with ρ=0
    W2 = pCN!(_W, 0.0)
    W2a = pCN(W, 0.0)
    WWrapper = NoiseWrapper(_W)

    # test source
    @test W2.source.W != W.W
    @test W2.source.Z == W.Z
    @test W2.source.W == W2a.source.W
    @test W2.source.Z == W2a.source.Z

    # statistics test
    W = WienerProcess(0.0, 0.0, 0.0)
    prob = NoiseProblem(W, (0.0, 1.0))
    sol = solve(prob, dt = 0.1)

    function prob_func(prob, i, repeat)
        _sol = deepcopy(sol)
        Wtmp = pCN!(_sol, 1.0)
        remake(prob, noise = Wtmp)
    end

    ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
    @time sim = solve(ensemble_prob, dt = 0.1, trajectories = 100)

    # Spot check the mean and the variance
    qs = 0:0.1:1
    for i in 1:11
        q = qs[i]
        @test ≈(timestep_mean(sim, i), sol.W[i], atol = 1e-2)
        @test ≈(timestep_meanvar(sim, i)[2], zero(sol.W[i]), atol = 1e-2)
    end

    # ρ = 0.5 test
    ρ = 0.5
    W = WienerProcess(0.0, 0.0, 0.0)
    prob = NoiseProblem(W, (0.0, 1000.0))
    sol = solve(prob, dt = 0.001)

    prob2 = NoiseProblem(pCN(sol, ρ), (0.0, 1000.0))
    sol2 = solve(prob2, dt = 0.001)

    function computedW(sim, indx)
        (sim.u[indx + 1] - sim.u[indx]) / sqrt(sim.t[indx + 1] - sim.t[indx])
    end
    dWnew = []
    dWold = []
    for i in 1:length(sol.t[1:(end - 1)])
        push!(dWnew, computedW(sol, i))
        push!(dWold, computedW(sol2, i))
    end
    @show cor(dWnew, dWold)

    @test ≈(cor(dWnew, dWold), ρ, atol = 1e-2)

    # Noise Grid tests

    trange = 0.0:0.001:1000.0

    # scalar
    Ws = cumsum([0.0;
        [sqrt(trange[i + 1] - ti) * randn()
         for (i, ti) in enumerate(trange[1:(end - 1)])]])
    NG = NoiseGrid(trange, Ws)
    NG2 = pCN(NG, 1.0) # ρ = 1.0
    NG3 = pCN(NG, 0.5) # ρ = 0.5

    @test ≈(cor(diff(NG.W) ./ sqrt.(diff(NG.t)), diff(NG2.W) ./ sqrt.(diff(NG2.t))), 1,
        atol = 1e-2)
    @test ≈(cor(diff(NG.W) ./ sqrt.(diff(NG.t)), diff(NG3.W) ./ sqrt.(diff(NG3.t))), 0.5,
        atol = 1e-2)

    # array
    dim = 4
    Ws = cumsum([[zeros(dim)];
        [sqrt(trange[i + 1] - ti) * randn(dim)
         for (i, ti) in enumerate(trange[1:(end - 1)])]])
    NG = NoiseGrid(trange, Ws)
    NG1 = pCN(NG, 1.0) # ρ = 1.0
    NG2 = pCN(NG, 0.5) # ρ = 0.5

    for i in 1:dim
        @test ≈(cor(getindex.(diff(NG.W) ./ sqrt.(diff(NG.t)), i),
                getindex.(diff(NG1.W) ./ sqrt.(diff(NG1.t)), i)), 1.0, atol = 1e-2)
        @test ≈(cor(getindex.(diff(NG.W) ./ sqrt.(diff(NG.t)), i),
                getindex.(diff(NG2.W) ./ sqrt.(diff(NG2.t)), i)), 0.5, atol = 1e-2)
    end
end
