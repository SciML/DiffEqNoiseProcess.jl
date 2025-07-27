using DiffEqNoiseProcess, Test, Random, RandomNumbers
using StochasticDiffEq, LinearAlgebra
@testset "NoiseWrapper" begin
    _W = WienerProcess(0.0, 0.0, 0.0)

    dt = 0.1
    calculate_step!(_W, dt, nothing, nothing)

    for i in 1:10
        accept_step!(_W, dt, nothing, nothing)
    end

    W2 = NoiseWrapper(_W)

    dt = 0.1
    calculate_step!(W2, dt, nothing, nothing)

    for i in 1:10
        accept_step!(W2, dt, nothing, nothing)
    end

    _W = WienerProcess(0.0, 0.0, 0.0)

    dt = 0.1
    calculate_step!(_W, dt, nothing, nothing)

    for i in 1:10
        accept_step!(_W, dt, nothing, nothing)
    end

    old_W = copy(_W.W)

    W2 = NoiseWrapper(_W)
    lspace = range(_W.t[1], stop = _W.t[end], length = 1000)
    dt = lspace[2] - lspace[1]
    calculate_step!(W2, dt, nothing, nothing)
    for t in lspace
        accept_step!(W2, dt, nothing, nothing)
    end

    @test W2.W[end] ≈ _W(W2.t[end])[1]
    @test W2.Z[end] ≈ _W(W2.t[end])[2]

    # Inplace

    _W = WienerProcess!(0.0, zeros(4), zeros(4))

    dt = 0.1
    calculate_step!(_W, dt, nothing, nothing)

    for i in 1:10
        accept_step!(_W, dt, nothing, nothing)
    end

    W2 = NoiseWrapper(_W)

    dt = 0.1
    calculate_step!(W2, dt, nothing, nothing)

    for i in 1:10
        accept_step!(W2, dt, nothing, nothing)
    end

    _W = WienerProcess!(0.0, zeros(4), zeros(4))

    dt = 0.1
    calculate_step!(_W, dt, nothing, nothing)

    for i in 1:10
        accept_step!(_W, dt, nothing, nothing)
    end

    W2 = NoiseWrapper(_W)
    lspace = range(_W.t[1], stop = _W.t[end], length = 20)
    dt = lspace[2] - lspace[1]
    calculate_step!(W2, dt, nothing, nothing)
    for t in lspace
        accept_step!(W2, dt, nothing, nothing)
    end

    @test W2.W[end] ≈ _W(W2.t[end])[1]
    @test W2.Z[end] ≈ _W(W2.t[end])[2]

    @test W2.W[end] ≈ _W.W[end - 1]
end

@testset "NoiseWrapper restart tests to interpolate" begin
    seed = 100
    u₀ = [0.75, 0.5]
    p = [-1.5, 0.05, 0.2, 0.01]
    trange = (0.0, 0.1)
    dtmix = trange[2] / 1e3

    function f_mixing!(du, u, p, t)
        du[1] = p[1] * u[1] + p[2] * u[2]
        du[2] = p[2] * u[1] + p[1] * u[2]
        nothing
    end

    function g_mixing!(du, u, p, t)
        du[1] = p[3] * u[1] + p[4] * u[2]
        du[2] = p[3] * u[1] + p[4] * u[2]
        nothing
    end

    function f_mixing(u, p, t)
        dx = p[1] * u[1] + p[2] * u[2]
        dy = p[2] * u[1] + p[1] * u[2]
        [dx, dy]
    end

    function g_mixing(u, p, t)
        dx = p[3] * u[1] + p[4] * u[2]
        dy = p[3] * u[1] + p[4] * u[2]
        [dx, dy]
    end

    Random.seed!(seed)
    prob = SDEProblem(f_mixing!, g_mixing!, u₀, trange, p)

    soltsave = collect(trange[1]:dtmix:trange[2])
    sol = solve(prob, EulerHeun(), dt = dtmix, save_noise = true, saveat = soltsave)

    Random.seed!(seed)
    proboop = SDEProblem(f_mixing, g_mixing, u₀, trange, p)
    soloop = solve(proboop, EulerHeun(), dt = dtmix, save_noise = true, saveat = soltsave)

    @test soloop.u ≈ sol.u

    interval = (sol.t[end - 1], sol.t[end])

    # oop

    _sol = deepcopy(soloop)
    soloop.W.save_everystep = false
    _sol.W.save_everystep = false

    forwardnoise = DiffEqNoiseProcess.NoiseWrapper(_sol.W, indx = 1000)
    checkWrapper = solve(
        remake(_sol.prob, tspan = interval, u0 = _sol(interval[1]),
            noise = forwardnoise),
        _sol.alg, save_noise = false;
        dt = abs(_sol.W.dt))

    @test checkWrapper.u[end - 1]≈soloop.u[end] rtol=1e-10
    @test checkWrapper.W.W[end]≈soloop.W.W[end] rtol=1e-16

    @show checkWrapper.u[end] - soloop.u[end]

    forwardnoise = DiffEqNoiseProcess.NoiseGrid(_sol.W.t[1000:1001], _sol.W.W[1000:1001])
    checkGrid = solve(
        remake(_sol.prob, tspan = interval, u0 = _sol(interval[1]),
            noise = forwardnoise),
        _sol.alg, save_noise = false;
        dt = abs(_sol.W.dt))

    @test checkGrid.u[end - 1]≈soloop.u[end] rtol=1e-10
    @test checkGrid.W.W[end]≈soloop.W.W[end] rtol=1e-16

    @show checkGrid.u[end] - soloop.u[end]

    # inplace

    _sol = deepcopy(sol)
    sol.W.save_everystep = false
    _sol.W.save_everystep = false

    forwardnoise = DiffEqNoiseProcess.NoiseWrapper(_sol.W, indx = 1000)
    checkWrapper = solve(
        remake(_sol.prob, tspan = interval, u0 = _sol(interval[1]),
            noise = forwardnoise),
        _sol.alg, save_noise = false;
        dt = abs(_sol.W.dt))

    @test checkWrapper.u[end - 1]≈sol.u[end] rtol=1e-10
    @test checkWrapper.W.W[end]≈sol.W.W[end] rtol=1e-16

    @show checkWrapper.u[end] - sol.u[end]

    forwardnoise = DiffEqNoiseProcess.NoiseGrid(_sol.W.t[1000:1001], _sol.W.W[1000:1001])
    checkGrid = solve(
        remake(_sol.prob, tspan = interval, u0 = _sol(interval[1]),
            noise = forwardnoise),
        _sol.alg, save_noise = false;
        dt = abs(_sol.W.dt))

    @test checkGrid.u[end - 1]≈sol.u[end] rtol=1e-10
    @test checkGrid.W.W[end]≈sol.W.W[end] rtol=1e-16

    @show checkGrid.u[end] - sol.u[end]
end

@testset "Interpolation with small eps jumps" begin
    seed = 100
    u₀ = [1.0]
    p = [1.1, 0.87]
    trange = (0.0, 1.0)
    dtmix = trange[2] / 1e4

    function f!(du, u, p, t)
        du[1] = p[1] * u[1]
        nothing
    end

    function g!(du, u, p, t)
        du[1] = p[2] * u[1]
        nothing
    end

    function f(u, p, t)
        dx = p[1] * u[1]
        [dx]
    end

    function g(u, p, t)
        dx = p[2] * u[1]
        [dx]
    end

    Random.seed!(seed)
    prob = SDEProblem(f!, g!, u₀, trange, p)
    sol = solve(prob, EM(), dt = dtmix, save_noise = true)

    Random.seed!(seed)
    proboop = SDEProblem(f, g, u₀, trange, p)
    soloop = solve(proboop, EM(), dt = dtmix, save_noise = true, save_end = false)

    @test soloop.u ≈ sol.u

    @test length(sol.u) != 1e4 + 1

    indx1 = Int(length(sol.t) - 100)
    interval = (sol.t[indx1], sol.t[end])

    #inplace

    _sol = deepcopy(sol)
    sol.W.save_everystep = false
    _sol.W.save_everystep = false

    forwardnoise = DiffEqNoiseProcess.NoiseWrapper(_sol.W, indx = indx1)

    checkWrapper = solve(
        remake(_sol.prob, tspan = interval, u0 = _sol(interval[1]),
            noise = forwardnoise),
        _sol.alg, save_noise = false;
        dt = abs(_sol.W.t[end - 1] - _sol.W.t[end - 2]),
        tstops = sol.t[indx1:end])

    @test length(checkWrapper.u) != length(sol.u)

    @test checkWrapper.u[1]≈sol.u[indx1] rtol=1e-10
    @test checkWrapper.u[end - 1]≈sol.u[end - 1] rtol=1e-10
    @test checkWrapper.u[end]≈sol.u[end] rtol=1e-10
    @test checkWrapper.W.W[end]≈sol.W.W[end] rtol=1e-10
    @test checkWrapper(sol.t[indx1:end]).u≈sol.u[indx1:end] rtol=1e-10

    _sol = deepcopy(soloop)
    sol.W.save_everystep = false
    _sol.W.save_everystep = false

    forwardnoise = DiffEqNoiseProcess.NoiseWrapper(_sol.W, indx = indx1)

    checkWrapper = solve(
        remake(_sol.prob, tspan = interval, u0 = _sol(interval[1]),
            noise = forwardnoise),
        _sol.alg, save_noise = false;
        dt = abs(_sol.W.t[end - 1] - _sol.W.t[end - 2]),
        tstops = sol.t[indx1:end])

    @test length(checkWrapper.u) != length(soloop.u)

    @test checkWrapper.u[1]≈sol.u[indx1] rtol=1e-10
    @test checkWrapper.u[end - 1]≈soloop.u[end - 1] rtol=1e-10
    @test checkWrapper.u[end]≈soloop.u[end] rtol=1e-10
    @test checkWrapper.W.W[end]≈soloop.W.W[end] rtol=1e-10
    @test checkWrapper(soloop.t[indx1:end]).u≈soloop.u[indx1:end] rtol=1e-10
end

@testset "Diagonal to Non-diagonal NoiseProcess conversion" begin
    t = 0.0
    dtleft = 0.01
    rand_prototype = [0.0 0.0 0.0; 0.0 0.0 0.0]
    save_noise = true
    _seed = 0x952197dfddfdce1f

    W = WienerProcess(t, rand_prototype,
        save_everystep = save_noise,
        rng = Xorshifts.Xoroshiro128Plus(_seed))
    prob = NoiseProblem(W, (0.0, 1.0))
    W2 = solve(prob; dt = 0.1)
    bW = DiffEqNoiseProcess.vec_NoiseProcess(W2)
    @test bW.dW isa AbstractVector && W.dW isa AbstractMatrix
end
