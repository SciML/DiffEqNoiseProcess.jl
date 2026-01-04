"""
    OrnsteinUhlenbeck{T1, T2, T3}

Parameters for the Ornstein-Uhlenbeck process.

# Fields
- `Θ`: Mean reversion rate (higher values mean faster reversion)
- `μ`: Long-term mean (the value the process reverts to)
- `σ`: Volatility/diffusion coefficient

The process follows the SDE: dX_t = Θ(μ - X_t)dt + σ dW_t
"""
struct OrnsteinUhlenbeck{T1, T2, T3}
    Θ::T1
    μ::T2
    σ::T3
end
"""
    (X::OrnsteinUhlenbeck)(dW, W, dt, u, p, t, rng)

Generate Ornstein-Uhlenbeck process increments using the exact distribution.

This uses the analytical solution for the OU process, providing exact samples
from the distribution rather than a numerical approximation.

# Arguments
- `dW`: Noise increment container
- `W`: Current noise process state
- `dt`: Time step
- `u`, `p`, `t`: State, parameters, and time (for compatibility)
- `rng`: Random number generator

# Returns
The increment dW for the OU process over time dt

# Reference
http://www.math.ku.dk/~susanne/StatDiff/Overheads1b.pdf
"""
function (X::OrnsteinUhlenbeck)(dW, W, dt, u, p, t, rng) #dist
    if dW isa AbstractArray
        rand_val = wiener_randn(rng, dW)
    else
        rand_val = wiener_randn(rng, typeof(dW))
    end
    drift = X.μ .+ (W.curW .- X.μ) .* exp.(-X.Θ * dt)
    diffusion = X.σ .* sqrt.(-expm1.(-2X.Θ * dt) ./ (2X.Θ))
    return drift .+ rand_val .* diffusion .- W.curW
end

#=
http://www.tandfonline.com/doi/pdf/10.1080/14697688.2014.941913?needAccess=true
http://www.tandfonline.com/doi/full/10.1080/14697688.2014.941913?src=recsys
=#

#=
(qX + r) = Θ(μ-x) = Θμ - Θx
q = -Θ
r = Θμ
https://arxiv.org/pdf/1011.0067.pdf page 18
note that in the paper there is a typo in the formula for σ^2
=#
"""
    ou_bridge(dW, ou, W, W0, Wh, q, h, u, p, t, rng)

Generate Ornstein-Uhlenbeck bridge interpolation between two points.

Provides exact sampling from an OU process conditioned on both endpoints,
useful for adaptive time-stepping and interpolation.

# Arguments
- `dW`: Noise increment container
- `ou`: OrnsteinUhlenbeck parameters
- `W`: Current noise process state
- `W0`: Starting value
- `Wh`: Target value at end of interval
- `q`: Interpolation parameter (0 to 1)
- `h`: Total time interval
- `u`, `p`, `t`: State, parameters, and time (for compatibility)
- `rng`: Random number generator

# Returns
Interpolated OU process value that maintains the correct distribution

# References
- http://www.tandfonline.com/doi/pdf/10.1080/14697688.2014.941913
- https://arxiv.org/pdf/1011.0067.pdf (page 18)
"""
function ou_bridge(dW, ou, W, W0, Wh, q, h, u, p, t, rng)
    return if dW isa AbstractArray
        rand_vec = wiener_randn(rng, dW)
        var = @. ou.σ^2 * sinh(ou.Θ * (h * (1.0 - q))) * (sinh(ou.Θ * (q * h))) /
            (ou.Θ * sinh(ou.Θ * h))
        @. (W0 - ou.μ) * (sinh(ou.Θ * (h * (1.0 - q))) / sinh(ou.Θ * h) - 1.0) +
            (Wh + W0 - ou.μ) * sinh(ou.Θ * q * h) / sinh(ou.Θ * h) + sqrt(var) * rand_vec

    else
        var = ou.σ^2 * sinh(ou.Θ * (h * (1.0 - q))) * (sinh(ou.Θ * (q * h))) /
            (ou.Θ * sinh(ou.Θ * h))
        (W0 - ou.μ) * (sinh(ou.Θ * (h * (1.0 - q))) / sinh(ou.Θ * h) - 1.0) +
            (Wh + W0 - ou.μ) * sinh(ou.Θ * q * h) / sinh(ou.Θ * h) +
            sqrt(var) * wiener_randn(rng, typeof(dW))
    end
end

"""
    ou_bridge!(rand_vec, ou, W, W0, Wh, q, h, u, p, t, rng)

Generate Ornstein-Uhlenbeck bridge interpolation in-place.

This is the in-place version of `ou_bridge`, modifying the provided array
rather than allocating a new one.

# Arguments
- `rand_vec`: Array to fill with interpolated values
- `ou`: OrnsteinUhlenbeck parameters
- `W`: Current noise process state
- `W0`: Starting value
- `Wh`: Target value at end of interval
- `q`: Interpolation parameter (0 to 1)
- `h`: Total time interval
- `u`, `p`, `t`: State, parameters, and time (for compatibility)
- `rng`: Random number generator

# Effects
Modifies rand_vec to contain the interpolated OU process values
"""
function ou_bridge!(rand_vec, ou, W, W0, Wh, q, h, u, p, t, rng)
    wiener_randn!(rng, rand_vec)
    return @.. rand_vec = (W0 - ou.μ) * (sinh(ou.Θ * (h * (1.0 - q))) / sinh(ou.Θ * h) - 1.0) +
        (Wh + W0 - ou.μ) * sinh(ou.Θ * q * h) / sinh(ou.Θ * h) +
        sqrt(
        ou.σ^2 * sinh(ou.Θ * (h * (1.0 - q))) * (sinh(ou.Θ * (q * h))) /
            (ou.Θ * sinh(ou.Θ * h))
    ) * rand_vec
end

@doc doc"""
a `Ornstein-Uhlenbeck` process, which is a Wiener process defined
by the stochastic differential equation

```math
dX_t = \theta (\mu - X_t) dt + \sigma dW_t
```

The `OrnsteinUhlenbeckProcess` is distribution exact (meaning, not a numerical
solution of the stochastic differential equation, but instead follows the exact
distribution properties). The constructor is:

```julia
OrnsteinUhlenbeckProcess(Θ,μ,σ,t0,W0,Z0=nothing;kwargs...)
OrnsteinUhlenbeckProcess!(Θ,μ,σ,t0,W0,Z0=nothing;kwargs...)
```
"""
function OrnsteinUhlenbeckProcess(Θ, μ, σ, t0, W0, Z0 = nothing; kwargs...)
    ou = OrnsteinUhlenbeck(Θ, μ, σ)
    return NoiseProcess{false}(
        t0, W0, Z0, ou,
        (
            rand_vec, W, W0, Wh, q, h, u, p, t,
            rng,
        ) -> ou_bridge(
            rand_vec,
            ou,
            W,
            W0,
            Wh,
            q,
            h,
            u,
            p,
            t,
            rng
        ); kwargs...
    )
end

"""
    OrnsteinUhlenbeck!{T1, T2, T3}

In-place version of OrnsteinUhlenbeck parameters.

# Fields
- `Θ`: Mean reversion rate (higher values mean faster reversion)
- `μ`: Long-term mean (the value the process reverts to)
- `σ`: Volatility/diffusion coefficient

The process follows the SDE: dX_t = Θ(μ - X_t)dt + σ dW_t
"""
struct OrnsteinUhlenbeck!{T1, T2, T3}
    Θ::T1
    μ::T2
    σ::T3
end

"""
    (X::OrnsteinUhlenbeck!)(rand_vec, W, dt, u, p, t, rng)

Generate Ornstein-Uhlenbeck process increments in-place using the exact distribution.

This is the in-place version that modifies the provided array rather than
allocating a new one.

# Arguments
- `rand_vec`: Array to fill with noise increments
- `W`: Current noise process state
- `dt`: Time step
- `u`, `p`, `t`: State, parameters, and time (for compatibility)
- `rng`: Random number generator

# Effects
Modifies rand_vec to contain the OU process increments
"""
function (X::OrnsteinUhlenbeck!)(rand_vec, W, dt, u, p, t, rng) #dist!
    wiener_randn!(rng, rand_vec)
    return @.. rand_vec = X.μ + (W.curW - X.μ) * exp(-X.Θ * dt) +
        rand_vec * X.σ * sqrt((-expm1.(-2 * X.Θ .* dt)) / (2 * X.Θ)) - W.curW
end

@doc doc"""
A `Ornstein-Uhlenbeck` process, which is a Wiener process defined
by the stochastic differential equation

```math
dX_t = \theta (\mu - X_t) dt + \sigma dW_t
```

The `OrnsteinUhlenbeckProcess` is distribution exact (meaning, not a numerical
solution of the stochastic differential equation, but instead follows the exact
distribution properties). The constructor is:

```julia
OrnsteinUhlenbeckProcess(Θ,μ,σ,t0,W0,Z0=nothing;kwargs...)
OrnsteinUhlenbeckProcess!(Θ,μ,σ,t0,W0,Z0=nothing;kwargs...)
```
"""
function OrnsteinUhlenbeckProcess!(Θ, μ, σ, t0, W0, Z0 = nothing; kwargs...)
    ou = OrnsteinUhlenbeck!(Θ, μ, σ)
    return NoiseProcess{true}(
        t0,
        W0,
        Z0,
        ou,
        (
            rand_vec, W, W0, Wh, q, h, u, p, t,
            rng,
        ) -> ou_bridge!(
            rand_vec,
            ou,
            W,
            W0,
            Wh,
            q,
            h,
            u,
            p,
            t,
            rng
        );
        kwargs...
    )
end
