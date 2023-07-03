using StaticArrays
@testset "Correlated Wiener Process" begin
    using DiffEqNoiseProcess, DiffEqBase, Test, Statistics

    # oop
    ρ = 0.3
    Γ = [1 ρ; ρ 1]
    W = CorrelatedWienerProcess(Γ, 0.0, zeros(2), zeros(2))

    dt = 0.1
    calculate_step!(W, dt, nothing, nothing)
    for i in 1:10
        accept_step!(W, dt, nothing, nothing)
    end

    prob = NoiseProblem(W, (0.0, 1.0))
    sol = solve(prob; dt = 0.1)

    # inplace
    W = CorrelatedWienerProcess!(Γ, 0.0, zeros(2), zeros(2))

    dt = 0.1
    calculate_step!(W, dt, nothing, nothing)
    for i in 1:10
        accept_step!(W, dt, nothing, nothing)
    end

    prob = NoiseProblem(W, (0.0, 1.0))
    sol = solve(prob; dt = 0.01)

    output_func = (sol, i) -> (sol.dW, false)
    ensemble_prob = EnsembleProblem(prob, output_func = output_func)
    @time sol = Array(solve(ensemble_prob, dt = dt, trajectories = 1_000_000))

    @test zero(mean(sol, dims = 2)[:])≈mean(sol, dims = 2)[:] atol=1e-2
    @test Γ≈cov(sol, dims = 2) / dt rtol=1e-2

    # with StaticArrays
    Γ = @SMatrix [1.0 ρ
        ρ 1.0]
    W = CorrelatedWienerProcess(Γ, 0.0, @SVector(zeros(2)), @SVector(zeros(2)))

    dt = 0.1
    calculate_step!(W, dt, nothing, nothing)
    for i in 1:10
        accept_step!(W, dt, nothing, nothing)
    end

    prob = NoiseProblem(W, (0.0, 1.0))
    sol = solve(prob; dt = 0.1)

    @test typeof(sol[1]) <: SArray
    @test typeof(sol[end]) <: SArray
end
