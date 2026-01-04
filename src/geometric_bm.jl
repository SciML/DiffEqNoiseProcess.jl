"""
    GeometricBrownianMotion{T1, T2}

Parameters for the geometric Brownian motion process.

# Fields
- `μ`: Drift parameter (expected return rate)
- `σ`: Volatility parameter (standard deviation of returns)

The process follows the SDE: dX_t = μX_t dt + σX_t dW_t

This is commonly used in financial models, particularly the Black-Scholes model.
"""
struct GeometricBrownianMotion{T1, T2}
    μ::T1
    σ::T2
end

"""
    (X::GeometricBrownianMotion)(dW, W, dt, u, p, t, rng)

Generate geometric Brownian motion increments using the exact distribution.

This uses the analytical solution for GBM, providing exact samples from
the log-normal distribution rather than a numerical approximation.

# Arguments
- `dW`: Noise increment container
- `W`: Current noise process state
- `dt`: Time step
- `u`, `p`, `t`: State, parameters, and time (for compatibility)
- `rng`: Random number generator

# Returns
The increment dW for the GBM process over time dt
"""
function (X::GeometricBrownianMotion)(dW, W, dt, u, p, t, rng) #dist
    drift = X.μ - (1 / 2) * X.σ^2
    if dW isa AbstractArray
        rand_val = wiener_randn(rng, dW)
    else
        rand_val = wiener_randn(rng, typeof(dW))
    end
    new_val = @.. exp(drift * dt + X.σ * sqrt(dt) * rand_val)
    return W.W[end] * (new_val - 1)
end

#=
u = q*h
t = h
s = 0

t-u = h-qh = (1-q)h
t-s = h
u - s = q*h
https://math.stackexchange.com/questions/412470/conditional-distribution-in-brownian-motion

=#
"""
    gbm_bridge(dW, gbm, W, W0, Wh, q, h, u, p, t, rng)

Generate geometric Brownian motion bridge interpolation between two points.

Provides exact sampling from a GBM process conditioned on both endpoints,
useful for adaptive time-stepping and interpolation.

# Arguments
- `dW`: Noise increment container
- `gbm`: GeometricBrownianMotion parameters
- `W`: Current noise process state
- `W0`: Starting value
- `Wh`: Target value at end of interval
- `q`: Interpolation parameter (0 to 1)
- `h`: Total time interval
- `u`, `p`, `t`: State, parameters, and time (for compatibility)
- `rng`: Random number generator

# Returns
Interpolated GBM process value that maintains the correct log-normal distribution

# Reference
https://math.stackexchange.com/questions/412470/conditional-distribution-in-brownian-motion
"""
function gbm_bridge(dW, gbm, W, W0, Wh, q, h, u, p, t, rng)
    if dW isa AbstractArray
        gbb_mean = @.. log(W0) + q * (log(W0 + Wh) - log(W0))
        gbb_std = sqrt((1 - q) * q * abs(h)) * gbm.σ

        x = gbb_mean + gbb_std * wiener_randn(rng, dW)
        return exp.(x) - W0
    else
        gbb_mean = log(W0) + q * (log(W0 + Wh) - log(W0))
        gbb_std = sqrt((1 - q) * q * abs(h)) * gbm.σ

        x = gbb_mean + gbb_std * wiener_randn(rng, typeof(dW))
        return exp(x) - W0
    end
end
"""
    gbm_bridge!(rand_vec, gbm, W, W0, Wh, q, h, u, p, t, rng)

Generate geometric Brownian motion bridge interpolation in-place.

This is the in-place version of `gbm_bridge`, modifying the provided array
rather than allocating a new one.

# Arguments
- `rand_vec`: Array to fill with interpolated values
- `gbm`: GeometricBrownianMotion parameters
- `W`: Current noise process state
- `W0`: Starting value
- `Wh`: Target value at end of interval
- `q`: Interpolation parameter (0 to 1)
- `h`: Total time interval
- `u`, `p`, `t`: State, parameters, and time (for compatibility)
- `rng`: Random number generator

# Effects
Modifies rand_vec to contain the interpolated GBM process values
"""
function gbm_bridge!(rand_vec, gbm, W, W0, Wh, q, h, u, p, t, rng)
    wiener_randn!(rng, rand_vec)
    @.. rand_vec = gbm.σ * sqrt((1 - q) * q * abs(h)) * rand_vec +
        (log(W0) + q * (log(W0 + Wh) - log(W0)))
    return @.. rand_vec = exp(rand_vec) - W0
end

@doc doc"""
A `GeometricBrownianMotion` process is a Wiener process with
constant drift `μ` and constant diffusion `σ`. I.e. this is the solution of the
stochastic differential equation

```math
dX_t = \mu X_t dt + \sigma X_t dW_t
```

The `GeometricBrownianMotionProcess` is distribution exact (meaning, not a numerical
solution of the stochastic differential equation, but instead follows the exact
distribution properties). It can be back interpolated exactly as well. The constructor is:

```julia
GeometricBrownianMotionProcess(μ,σ,t0,W0,Z0=nothing;kwargs...)
GeometricBrownianMotionProcess!(μ,σ,t0,W0,Z0=nothing;kwargs...)
```
"""
function GeometricBrownianMotionProcess(μ, σ, t0, W0, Z0 = nothing; kwargs...)
    gbm = GeometricBrownianMotion(μ, σ)
    return NoiseProcess{false}(
        t0, W0, Z0, gbm,
        (
            dW, W, W0, Wh, q, h, u, p, t,
            rng,
        ) -> gbm_bridge(
            dW, gbm, W, W0,
            Wh, q, h, u, p, t,
            rng
        ); kwargs...
    )
end

"""
    GeometricBrownianMotion!{T1, T2}

In-place version of GeometricBrownianMotion parameters.

# Fields
- `μ`: Drift parameter (expected return rate)
- `σ`: Volatility parameter (standard deviation of returns)

The process follows the SDE: dX_t = μX_t dt + σX_t dW_t
"""
struct GeometricBrownianMotion!{T1, T2}
    μ::T1
    σ::T2
end
"""
    (X::GeometricBrownianMotion!)(rand_vec, W, dt, u, p, t, rng)

Generate geometric Brownian motion increments in-place using the exact distribution.

This is the in-place version that modifies the provided array rather than
allocating a new one.

# Arguments
- `rand_vec`: Array to fill with noise increments
- `W`: Current noise process state
- `dt`: Time step
- `u`, `p`, `t`: State, parameters, and time (for compatibility)
- `rng`: Random number generator

# Effects
Modifies rand_vec to contain the GBM process increments
"""
function (X::GeometricBrownianMotion!)(rand_vec, W, dt, u, p, t, rng) #dist!
    wiener_randn!(rng, rand_vec)
    return @.. rand_vec = W.W[end] * expm1(X.μ - (1 / 2) * X.σ * dt + X.σ * sqrt(dt) * rand_vec)
end

@doc doc"""
A `GeometricBrownianMotion` process is a Wiener process with
constant drift `μ` and constant diffusion `σ`. I.e. this is the solution of the
stochastic differential equation

```math
dX_t = \mu X_t dt + \sigma X_t dW_t
```

The `GeometricBrownianMotionProcess` is distribution exact (meaning, not a numerical
solution of the stochastic differential equation, but instead follows the exact
distribution properties). It can be back interpolated exactly as well. The constructor is:

```julia
GeometricBrownianMotionProcess(μ,σ,t0,W0,Z0=nothing;kwargs...)
GeometricBrownianMotionProcess!(μ,σ,t0,W0,Z0=nothing;kwargs...)
```
"""
function GeometricBrownianMotionProcess!(μ, σ, t0, W0, Z0 = nothing; kwargs...)
    gbm = GeometricBrownianMotion!(μ, σ)
    return NoiseProcess{true}(
        t0, W0, Z0, gbm,
        (
            rand_vec, W, W0, Wh, q, h, u, p, t,
            rng,
        ) -> gbm_bridge!(
            rand_vec,
            gbm, W, W0,
            Wh, q, h, u,
            p, t, rng
        );
        kwargs...
    )
end
