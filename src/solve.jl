"""
    solve(prob::AbstractNoiseProblem; dt, kwargs...)

Solve a noise problem by simulating the noise process over the specified time span.

## Arguments
- `prob`: An `AbstractNoiseProblem` containing the noise process and time span
- `dt`: The time step size (required, no default)

## Keyword Arguments
- `seed`: Random seed for reproducible results
- Additional keyword arguments are passed to the underlying solver

## Returns
A `NoiseProcess` object containing the simulated noise trajectory.

## Example
```julia
W = WienerProcess(0.0, 0.0, 1.0)
prob = NoiseProblem(W, (0.0, 1.0))
sol = solve(prob; dt = 0.01)
```
"""
function DiffEqBase.__solve(
        prob::AbstractNoiseProblem,
        args::Union{Nothing, SciMLBase.DEAlgorithm}...; dt = 0.0,
        kwargs...
    )
    if dt == 0.0 || dt == nothing
        error("dt must be provided to simulate a noise process. Please pass dt=...")
    end
    W = copy(prob.noise)
    if W isa Union{NoiseProcess, NoiseTransport}
        if prob.seed != 0
            Random.seed!(W.rng, prob.seed)
        else
            Random.seed!(W.rng, rand(UInt64))
        end
    end
    if W.reset
        reinit!(W, dt, t0 = prob.tspan[1])
    end

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
    return W
end
