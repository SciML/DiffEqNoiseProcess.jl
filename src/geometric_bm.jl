struct GeometricBrownianMotion{T1, T2}
    μ::T1
    σ::T2
end

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
function gbm_bridge!(rand_vec, gbm, W, W0, Wh, q, h, u, p, t, rng)
    wiener_randn!(rng, rand_vec)
    @.. rand_vec = gbm.σ * sqrt((1 - q) * q * abs(h)) * rand_vec +
                   (log(W0) + q * (log(W0 + Wh) - log(W0)))
    @.. rand_vec = exp(rand_vec) - W0
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
    NoiseProcess{false}(t0, W0, Z0, gbm,
        (dW, W, W0, Wh, q, h, u, p, t,
            rng) -> gbm_bridge(dW, gbm, W, W0,
            Wh, q, h, u, p, t,
            rng); kwargs...)
end

struct GeometricBrownianMotion!{T1, T2}
    μ::T1
    σ::T2
end
function (X::GeometricBrownianMotion!)(rand_vec, W, dt, u, p, t, rng) #dist!
    wiener_randn!(rng, rand_vec)
    @.. rand_vec = W.W[end] * expm1(X.μ - (1 / 2) * X.σ * dt + X.σ * sqrt(dt) * rand_vec)
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
    NoiseProcess{true}(t0, W0, Z0, gbm,
        (rand_vec, W, W0, Wh, q, h, u, p, t,
            rng) -> gbm_bridge!(rand_vec,
            gbm, W, W0,
            Wh, q, h, u,
            p, t, rng);
        kwargs...)
end
