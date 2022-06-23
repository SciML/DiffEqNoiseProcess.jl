function DiffEqBase.__solve(
    prob::AbstractNoiseProblem,
    args::Union{Nothing,SciMLBase.DEAlgorithm}...;
    dt = 0.0,
    kwargs...,
)
    if dt == 0.0 || dt == nothing
        error("dt must be provided to simulate a noise process. Please pass dt=...")
    end
    W = deepcopy(prob.noise)
    if typeof(W) <: NoiseProcess
        if prob.seed != 0
            Random.seed!(W.rng, prob.seed)
        else
            Random.seed!(W.rng, rand(UInt64))
        end
    end
    W.curt = prob.tspan[1]
    W.dt = dt
    setup_next_step!(W, nothing, nothing)
    tType = typeof(W.curt)
    while W.curt < prob.tspan[2]
        if tType <: AbstractFloat && abs(W.curt + dt - prob.tspan[2]) < 100 * eps(dt) # Correct the end due to floating point error
            dt = prob.tspan[2] - W.curt
            accept_step!(W, dt, nothing, nothing)
            W.curt = prob.tspan[2]
        else
            accept_step!(W, dt, nothing, nothing)
        end
    end
    W
end
