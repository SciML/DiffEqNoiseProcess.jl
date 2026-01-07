"""
    solve(prob::AbstractNoiseProblem; dt, saveat=nothing, kwargs...)

Solve a noise problem by simulating the noise process over the specified time span.

## Arguments
- `prob`: An `AbstractNoiseProblem` containing the noise process and time span
- `dt`: The time step size (required, no default)

## Keyword Arguments
- `seed`: Random seed for reproducible results
- `saveat`: A collection of time points at which to save the noise process values.
  If provided, the noise will be interpolated at these times and only those values
  will be stored. This is useful when an SDE uses `saveat` and you want to extract
  the noise process only at those times without allocating memory for the entire
  process history.
- Additional keyword arguments are passed to the underlying solver

## Returns
A `NoiseProcess` object containing the simulated noise trajectory.

## Example
```julia
W = WienerProcess(0.0, 0.0, 1.0)
prob = NoiseProblem(W, (0.0, 1.0))
sol = solve(prob; dt = 0.01)

# Save only at specific times
sol = solve(prob; dt = 0.01, saveat = 0.0:0.1:1.0)
```
"""
function DiffEqBase.__solve(
        prob::AbstractNoiseProblem,
        args::Union{Nothing, SciMLBase.DEAlgorithm}...; dt = 0.0,
        saveat = nothing,
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

    # When saveat is provided, we need to save all steps during simulation
    # so we can interpolate later, regardless of the save_everystep setting
    original_save_everystep = nothing
    if saveat !== nothing && hasproperty(W, :save_everystep)
        original_save_everystep = W.save_everystep
        W.save_everystep = true
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

    # Handle saveat: interpolate at specified times and replace stored values
    if saveat !== nothing
        _apply_saveat!(W, saveat)
    end

    return W
end

"""
    _apply_saveat!(W, saveat)

Internal function to apply saveat to a noise process after simulation.
Interpolates the noise process at the specified times and replaces the stored
values with only those at the saveat times.
"""
function _apply_saveat!(W, saveat)
    # Convert saveat to a sorted vector if needed
    saveat_vec = collect(saveat)

    # Store original save_everystep setting and temporarily disable it
    # to prevent interpolate! from inserting additional values
    original_save_everystep = nothing
    if hasproperty(W, :save_everystep)
        original_save_everystep = W.save_everystep
        W.save_everystep = false
    end

    # Collect interpolated values at each saveat time
    new_t = similar(W.t, 0)
    new_W = similar(W.W, 0)
    has_Z = W.Z !== nothing
    if has_Z
        new_Z = similar(W.Z, 0)
    end

    t_min = W.t[1]
    t_max = W.t[end]

    for t in saveat_vec
        # Skip times outside the simulation range
        if t < t_min || t > t_max
            continue
        end
        # Interpolate at time t
        W_val, Z_val = interpolate!(W, nothing, nothing, t)
        push!(new_t, t)
        push!(new_W, copy(W_val))
        if has_Z
            push!(new_Z, copy(Z_val))
        end
    end

    # Replace stored values with saveat values
    resize!(W.t, length(new_t))
    resize!(W.W, length(new_W))
    W.t .= new_t
    W.W .= new_W
    W.u = W.W  # Keep alias consistent

    if has_Z
        resize!(W.Z, length(new_Z))
        W.Z .= new_Z
    end

    # Restore save_everystep setting
    if original_save_everystep !== nothing
        W.save_everystep = original_save_everystep
    end

    return nothing
end
