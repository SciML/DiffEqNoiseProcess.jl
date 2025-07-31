struct OrnsteinUhlenbeck{T1, T2, T3}
    Θ::T1
    μ::T2
    σ::T3
end
# http://www.math.ku.dk/~susanne/StatDiff/Overheads1b.pdf
function (X::OrnsteinUhlenbeck)(dW, W, dt, u, p, t, rng) #dist
    if dW isa AbstractArray
        rand_val = wiener_randn(rng, dW)
    else
        rand_val = wiener_randn(rng, typeof(dW))
    end
    drift = X.μ .+ (W.curW .- X.μ) .* exp.(-X.Θ * dt)
    diffusion = X.σ .* sqrt.(-expm1.(-2X.Θ * dt) ./ (2X.Θ))
    drift .+ rand_val .* diffusion .- W.curW
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
function ou_bridge(dW, ou, W, W0, Wh, q, h, u, p, t, rng)
    if dW isa AbstractArray
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

function ou_bridge!(rand_vec, ou, W, W0, Wh, q, h, u, p, t, rng)
    wiener_randn!(rng, rand_vec)
    @.. rand_vec = (W0 - ou.μ) * (sinh(ou.Θ * (h * (1.0 - q))) / sinh(ou.Θ * h) - 1.0) +
                   (Wh + W0 - ou.μ) * sinh(ou.Θ * q * h) / sinh(ou.Θ * h) +
                   sqrt(ou.σ^2 * sinh(ou.Θ * (h * (1.0 - q))) * (sinh(ou.Θ * (q * h))) /
                        (ou.Θ * sinh(ou.Θ * h))) * rand_vec
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
    NoiseProcess{false}(t0, W0, Z0, ou,
        (rand_vec, W, W0, Wh, q, h, u, p, t,
            rng) -> ou_bridge(rand_vec,
            ou,
            W,
            W0,
            Wh,
            q,
            h,
            u,
            p,
            t,
            rng); kwargs...)
end

struct OrnsteinUhlenbeck!{T1, T2, T3}
    Θ::T1
    μ::T2
    σ::T3
end

function (X::OrnsteinUhlenbeck!)(rand_vec, W, dt, u, p, t, rng) #dist!
    wiener_randn!(rng, rand_vec)
    @.. rand_vec = X.μ + (W.curW - X.μ) * exp(-X.Θ * dt) +
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
    NoiseProcess{true}(t0,
        W0,
        Z0,
        ou,
        (rand_vec, W, W0, Wh, q, h, u, p, t,
            rng) -> ou_bridge!(rand_vec,
            ou,
            W,
            W0,
            Wh,
            q,
            h,
            u,
            p,
            t,
            rng);
        kwargs...)
end
