"""
    cpp_bridge(dW, cpp, W, W0, Wh, q, h, u, p, t, rng)

Generate compound Poisson process bridge interpolation between two points.

Uses binomial thinning to distribute jumps appropriately between endpoints.

# Arguments
- `dW`: Noise increment container
- `cpp`: CompoundPoissonProcess parameters
- `W`: Current noise process state
- `W0`: Starting value
- `Wh`: Jump count difference (must be integer)
- `q`: Interpolation parameter (0 to 1)
- `h`: Total time interval
- `u`, `p`, `t`: State, parameters, and time (for compatibility)
- `rng`: Random number generator

# Returns
Number of jumps distributed according to binomial thinning
"""
function cpp_bridge(dW, cpp, W, W0, Wh, q, h, u, p, t, rng)
    rand.(rng, Distributions.Binomial.(Int.(Wh), float.(q)))
end
"""
    cpp_bridge!(rand_vec, cpp, W, W0, Wh, q, h, u, p, t, rng)

In-place version of `cpp_bridge`.

# Effects
Modifies rand_vec to contain binomially distributed jump counts
"""
function cpp_bridge!(rand_vec, cpp, W, W0, Wh, q, h, u, p, t, rng)
    rand_vec .= rand.(rng, Distributions.Binomial.(Int.(Wh), float.(q)))
end

"""
    CompoundPoissonProcess{R, CR}

A compound Poisson process for modeling jump processes.

The process has jumps that occur according to a Poisson process with given rate,
and jump sizes determined by a specified distribution.

# Fields
- `rate`: Jump rate function or constant (λ parameter)
- `currate`: Current rate value (cached for efficiency)
- `computerates`: Whether to recompute rates at each step

# Constructor
```julia
CompoundPoissonProcess(rate, t0, W0; computerates = true, kwargs...)
```

# Examples
```julia
# Constant rate
proc = CompoundPoissonProcess(2.0, 0.0, 0.0)

# Time-dependent rate
rate_func(u, p, t) = 1.0 + 0.5*sin(t)
proc = CompoundPoissonProcess(rate_func, 0.0, 0.0)
```

# References
https://www.math.wisc.edu/~anderson/papers/AndPostleap.pdf
Incorporating postleap checks in tau-leaping
J. Chem. Phys. 128, 054103 (2008); https://doi.org/10.1063/1.2819665
"""
mutable struct CompoundPoissonProcess{R, CR}
    rate::R
    currate::CR
    computerates::Bool
    function CompoundPoissonProcess(rate, t0, W0; computerates = true, kwargs...)
        cpp = new{typeof(rate), typeof(W0)}(rate, W0, computerates)
        NoiseProcess{false}(t0, W0, nothing, cpp,
            (dW, W, W0, Wh, q, h, u, p, t, rng) -> cpp_bridge(dW, cpp, W,
                W0, Wh, q, h,
                u, p, t, rng);
            continuous = false, cache = cpp, kwargs...)
    end
end
"""
    (P::CompoundPoissonProcess)(dW, W, dt, u, p, t, rng)

Generate compound Poisson process increments.

Samples from a Poisson distribution with rate λ*dt to determine the number
of jumps in the time interval dt.

# Arguments
- `dW`: Noise increment container
- `W`: Current noise process state
- `dt`: Time step
- `u`: Current state (for state-dependent rates)
- `p`: Parameters
- `t`: Current time
- `rng`: Random number generator

# Returns
Number of jumps in the interval dt
"""
function (P::CompoundPoissonProcess)(dW, W, dt, u, p, t, rng)
    P.computerates && (P.currate = P.rate(u, p, t))
    PoissonRandom.pois_rand.(rng, dt .* P.currate)
end

"""
    CompoundPoissonProcess!{R, CR}

In-place version of CompoundPoissonProcess.

See `CompoundPoissonProcess` for details.

# Constructor
```julia
CompoundPoissonProcess!(rate, t0, W0; computerates = true, kwargs...)
```
"""
struct CompoundPoissonProcess!{R, CR}
    rate::R
    currate::CR
    computerates::Bool
    function CompoundPoissonProcess!(rate, t0, W0; computerates = true, kwargs...)
        cpp = new{typeof(rate), typeof(W0)}(rate, copy(W0), computerates)
        NoiseProcess{true}(t0, W0, nothing, cpp,
            (rand_vec, W, W0, Wh, q, h, u, p, t, rng) -> cpp_bridge!(rand_vec,
                cpp, W,
                W0, Wh,
                q, h, u,
                p, t,
                rng);
            continuous = false, cache = cpp, kwargs...)
    end
end
"""
    (P::CompoundPoissonProcess!)(rand_vec, W, dt, u, p, t, rng)

Generate compound Poisson process increments in-place.

This is the in-place version that modifies the provided array.

# Arguments
- `rand_vec`: Array to fill with jump counts
- `W`: Current noise process state
- `dt`: Time step
- `u`: Current state (for state-dependent rates)
- `p`: Parameters
- `t`: Current time
- `rng`: Random number generator

# Effects
Modifies rand_vec to contain Poisson-distributed jump counts
"""
function (P::CompoundPoissonProcess!)(rand_vec, W, dt, u, p, t, rng)
    P.computerates && P.rate(P.currate, u, p, t)
    @. rand_vec = PoissonRandom.pois_rand(rng, dt * P.currate)
end
