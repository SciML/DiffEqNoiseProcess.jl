@testset "SDE Adaptive Distribution Tests" begin
    using StochasticDiffEq, StatsBase, Distributions, HypothesisTests
    using Random
    using DiffEqProblemLibrary.SDEProblemLibrary
    # load problems
    SDEProblemLibrary.importsdeproblems()

    prob = prob_sde_linear
    Random.seed!(200)
    N = 100
    M = 5
    ps = Vector{Float64}(undef, M)
    T = prob.tspan[2]

    for j in 1:M
        Wends = Vector{Float64}(undef, N)
        for i in 1:N
            W = WienerProcess(0.0, 0.0, 0.0, rswm = RSWM(adaptivealg = :RSwM1))
            _prob = remake(prob, noise = W)
            sol = solve(_prob, SRI(), dt = 1 / 2^(4), abstol = 1e-2, reltol = 0)
            Wends[i] = sol.W.W[end]
        end
        kssol = ApproximateOneSampleKSTest(Wends / sqrt(T), Normal())
        ps[j] = pvalue(kssol) #Should be not significant (most of the time)
    end

    @test sum(ps .> 0.05) > length(ps) / 2 ### Make sure more passes than fails

    for j in 1:M
        Wends = Vector{Float64}(undef, N)
        for i in 1:N
            W = WienerProcess(0.0, 0.0, 0.0, rswm = RSWM(adaptivealg = :RSwM2))
            _prob = remake(prob, noise = W)
            sol = solve(_prob, SRI(), dt = 1 / 2^(4), abstol = 1e-2, reltol = 0)
            Wends[i] = sol.W.W[end]
        end
        kssol = ApproximateOneSampleKSTest(Wends / sqrt(T), Normal())
        ps[j] = pvalue(kssol) #Should be not significant (most of the time)
    end

    @test sum(ps .> 0.05) > length(ps) / 2 ### Make sure more passes than fails

    for j in 1:M
        Wends = Vector{Float64}(undef, N)
        for i in 1:N
            W = WienerProcess(0.0, 0.0, 0.0, rswm = RSWM(adaptivealg = :RSwM3))
            _prob = remake(prob, noise = W)
            sol = solve(_prob, SRI(), dt = 1 / 2^(4), abstol = 1e-2, reltol = 0)
            Wends[i] = sol.W.W[end]
        end
        kssol = ApproximateOneSampleKSTest(Wends / sqrt(T), Normal())
        ps[j] = pvalue(kssol) #Should be not significant (most of the time)
    end

    @test sum(ps .> 0.05) > length(ps) / 2 ### Make sure more passes than fails
end
