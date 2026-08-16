function construct_correlated_noisefunc(Γ)
    γ = svd(Γ)
    A = γ.U * Diagonal(sqrt.(γ.S))
    dist = function (dW, W, dt, u, p, t, rng)
        if dW isa AbstractArray
            return A * sqrt.(abs(dt)) * wiener_randn(rng, dW)
        else
            return A * sqrt.(abs(dt)) * wiener_randn(rng, typeof(dW))
        end
    end
    bridge = function (dW, W, W0, Wh, q, h, u, p, t, rng)
        error("Bridging distribution is unknown. Cannot use adapativity")
    end
    return dist, bridge
end

"""
    CorrelatedWienerProcess(Γ, t0, W0, Z0 = nothing; rng = Random.default_rng())

Construct an out-of-place Wiener process with constant covariance matrix `Γ`.
The covariance is factored once and applied to independent normal increments.

```julia
Γ = [1.0 0.2; 0.2 1.0]
W = CorrelatedWienerProcess(Γ, 0.0, zeros(2))
```

# Arguments
- `Γ`: Square, positive-semidefinite covariance matrix.
- `t0`: Initial time.
- `W0`: Initial process value and prototype; its length must match `Γ`.
- `Z0`: Optional auxiliary process value.

# Keywords
- `rng`: Random number generator used for increments.

# Returns
A `NoiseProcess` with `isinplace(W) == false`.
"""
function CorrelatedWienerProcess(
        Γ, t0, W0, Z0 = nothing;
        rng = Random.default_rng()
    )
    return NoiseProcess{false}(
        t0, W0, Z0, construct_correlated_noisefunc(Γ)..., rswm = RSWM(),
        rng = rng, covariance = Γ
    )
end

function construct_correlated_noisefunc!(Γ)
    γ = svd(Γ)
    A = γ.U * Diagonal(sqrt.(γ.S))
    # The scratch vector must live inside the call: the closure is shared across
    # `copy`s of the process, so a captured buffer is raced on by ensemble
    # trajectories running on different threads.
    dist = function (rand_vec, W, dt, u, p, t, rng)
        b = Vector{eltype(rand_vec)}(undef, length(rand_vec))
        wiener_randn!(rng, b)
        b .*= sqrt.(abs(dt))
        return mul!(rand_vec, A, b)
    end
    bridge = function (rand_vec, W, W0, Wh, q, h, u, p, t, rng)
        error("Bridging distribution is unknown. Cannot use adapativity")
    end
    return dist, bridge
end

"""
    CorrelatedWienerProcess!(Γ, t0, W0, Z0 = nothing; rng = Random.default_rng())

Construct the in-place variant of [`CorrelatedWienerProcess`](@ref). It uses
the same constant covariance model and writes increments into `W0`-shaped
storage.

```julia
Γ = [1.0 0.2; 0.2 1.0]
W = CorrelatedWienerProcess!(Γ, 0.0, zeros(2))
calculate_step!(W, 0.01, nothing, nothing)
```

# Arguments
- `Γ`: Square, positive-semidefinite covariance matrix.
- `t0`: Initial time.
- `W0`: Initial process value and storage prototype.
- `Z0`: Optional auxiliary process value.

# Keywords
- `rng`: Random number generator used for increments.

# Returns
A `NoiseProcess` with `isinplace(W) == true`.
"""
function CorrelatedWienerProcess!(
        Γ, t0, W0, Z0 = nothing;
        rng = Random.default_rng()
    )
    return NoiseProcess{true}(
        t0, W0, Z0, construct_correlated_noisefunc!(Γ)...,
        rswm = RSWM(), rng = rng, covariance = Γ
    )
end
