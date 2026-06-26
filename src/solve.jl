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
# After the endpoint is snapped onto `tspan[2]`, the next step prepared inside
# `accept_step!` still corresponds to the unsnapped `curt`. For function-evaluated
# noises the prepared step is a pure function of `curt`, so recompute it at the
# snapped time to keep `dW == W(curt + dt) - curW` exact. Stateful processes can't
# be re-prepared (RSwM stack pops, grid domain guards), and for them the prepared
# step is random so the sub-eps time label drift is unobservable.
resync_prepared_step!(W::Union{NoiseFunction, NoiseTransport}) =
    calculate_step!(W, W.dt, nothing, nothing)
resync_prepared_step!(W) = nothing

function DiffEqBase.__solve(
        prob::AbstractNoiseProblem,
        args::Union{Nothing, SciMLBase.AbstractDEAlgorithm}...; dt = 0.0,
        kwargs...
    )
    if dt == 0.0 || dt == nothing
        error("dt must be provided to simulate a noise process. Please pass dt=...")
    end
    W = copy(prob.noise)
    if W isa Union{NoiseProcess, NoiseTransport}
        if prob.seed != 0
            Random.seed!(W.rng, prob.seed)
        elseif W.reseed
            Random.seed!(W.rng)
        end
    end
    if W.reset
        reinit!(W, dt, t0 = prob.tspan[1])
    end

    setup_next_step!(W, nothing, nothing)
    tType = typeof(W.curt)
    dtType = typeof(dt)
    while W.curt < prob.tspan[2]
        # Tolerance scaled by the magnitude of the time values rather than by `dt`:
        # the floating point drift accumulated over many steps is on the order of
        # `eps(tspan[2])`, which can be far larger than `eps(dt)` when `dt` is small.
        # The previous `100 * eps(dt)` tolerance was therefore essentially never hit,
        # so the solve took a full extra step and stopped past `tspan[2]`.
        #
        # When `dt` is given at a coarser precision than the time type (e.g. a Float32
        # `dt = 1/50f0` with a Float64 tspan), the per-step error is set by `dt`'s
        # precision, so the drift accumulated over the span reaches ~`eps(dtType(tspan[2]))`
        # (~1e-7 for Float32) — orders of magnitude larger than `eps(tType(tspan[2]))`
        # (~1e-16 for Float64). Use the coarser of the two so the near-endpoint step is
        # snapped onto `tspan[2]` instead of having a sub-`dt` sliver step appended past it.
        mag = max(abs(prob.tspan[2]), abs(W.curt))
        endtol = tType <: AbstractFloat ?
            100 * max(
                eps(tType(mag)),
                dtType <: AbstractFloat ? eps(dtType(mag)) : zero(tType(mag))
            ) :
            zero(W.curt)
        if tType <: AbstractFloat && abs(prob.tspan[2] - (W.curt + W.dt)) <= endtol
            # The prepared step lands on `tspan[2]` up to floating point drift. Take it as
            # usual but snap the recorded endpoint exactly onto `tspan[2]` so the solution
            # neither overshoots nor stops just short of the requested final time.
            accept_step!(W, dt, nothing, nothing)
            W.curt = prob.tspan[2]
            # Only NoiseProcess/SimpleNoiseProcess record a `t` timeseries to snap;
            # types like NoiseFunction and NoiseGrid have no `save_everystep` field.
            hasfield(typeof(W), :save_everystep) && W.save_everystep &&
                (W.t[end] = prob.tspan[2])
            resync_prepared_step!(W)
        elseif tType <: AbstractFloat && W.curt + W.dt > prob.tspan[2] + endtol
            # The prepared step would overshoot `tspan[2]` by more than rounding. Recompute
            # it at exactly the remaining width so the solution ends on `tspan[2]`.
            dtcorrect = prob.tspan[2] - W.curt
            calculate_step!(W, dtcorrect, nothing, nothing)
            accept_step!(W, dtcorrect, nothing, nothing)
            W.curt = prob.tspan[2]
            hasfield(typeof(W), :save_everystep) && W.save_everystep &&
                (W.t[end] = prob.tspan[2])
            resync_prepared_step!(W)
        else
            accept_step!(W, dt, nothing, nothing)
        end
    end
    return W
end
