@testset "Remake" begin
    using SciMLBase, DiffEqNoiseProcess, Test, Random
    W = WienerProcess(0.0, 1.0, 1.0, rng = Xoshiro(42))
    dt = 0.1
    W.dt = dt
    u = nothing
    p = nothing # for state-dependent distributions
    calculate_step!(W, dt, u, p)
    for i in 1:10
        accept_step!(W, dt, u, p)
    end
    W2 = copy(W)
    for prop in propertynames(W)
        @test getfield(W, prop) === getfield(W2, prop)
    end
    rng2 = Xoshiro(43)
    W3 = remake(W2, rng = rng2)
    @test W3.rng === rng2
    W.rng = rng2
    for prop in propertynames(W)
        @test getfield(W, prop) === getfield(W3, prop)
    end
end
