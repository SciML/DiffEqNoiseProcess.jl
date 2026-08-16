"""
    BrownianBridge(t0, tend, W0, Wend, Z0 = nothing, Zend = nothing; kwargs...)

Construct a Wiener process conditioned on its values at `t0` and `tend`.
The process is distribution exact and can be interpolated between the endpoints.

```julia
bridge = BrownianBridge(0.0, 1.0, 0.0, 1.0)
bridge(0.5)
```

# Arguments
- `t0`, `tend`: Start and end times.
- `W0`, `Wend`: Process values at the two endpoints.
- `Z0`, `Zend`: Optional endpoint values for the auxiliary process.

# Keywords
Additional keywords are forwarded to [`WienerProcess`](@ref).

# Returns
A `NoiseProcess` initialized with the endpoint increment in its RSwM stack.
"""
function BrownianBridge(t0, tend, W0, Wend, Z0 = nothing, Zend = nothing; kwargs...)
    W = WienerProcess(t0, W0, Z0; kwargs...)
    h = tend - t0
    Wh = Wend - W0
    if Z0 !== nothing
        Zh = Zend - Z0
    else
        Zh = nothing
    end
    push!(W.S₁, (h, Wh, Zh))
    push!(W.reinitS₁, (h, Wh, Zh))
    return W
end

"""
    BrownianBridge!(t0, tend, W0, Wend, Z0 = nothing, Zend = nothing; kwargs...)

Construct the in-place variant of [`BrownianBridge`](@ref). The endpoint arrays
are modified when their increments are formed, so pass copies when the original
endpoint values must be preserved.

```julia
bridge = BrownianBridge!(0.0, 1.0, zeros(2), ones(2))
bridge(0.5)
```

# Arguments
- `t0`, `tend`: Start and end times.
- `W0`, `Wend`: Array-valued endpoint values.
- `Z0`, `Zend`: Optional auxiliary endpoint values.

# Keywords
Additional keywords are forwarded to [`WienerProcess!`](@ref).

# Returns
A `NoiseProcess` with `isinplace(W) == true`.
"""
function BrownianBridge!(t0, tend, W0, Wh, Z0 = nothing, Zh = nothing; kwargs...)
    W = WienerProcess!(t0, W0, Z0; kwargs...)
    h = tend - t0
    Wh .-= W0
    if Z0 !== nothing
        Zh .-= Z0
    else
        Zh = nothing
    end
    push!(W.S₁, (h, Wh, Zh))
    push!(W.reinitS₁, (h, Wh, Zh))
    return W
end

"""
    GeometricBrownianBridge(μ, σ, t0, tend, W0, Wend,
        Z0 = nothing, Zend = nothing; kwargs...)

A `GeometricBrownianBridge` is a geometric Brownian motion process with
pre-defined start and end values.

This creates a GBM process that is conditioned to pass through specific values at the beginning and end of the time interval,
useful for financial modeling where asset prices must match observed values.

# Arguments
- `μ`: Drift parameter
- `σ`: Volatility parameter
- `t0`: Starting time
- `tend`: Ending time
- `W0`: Starting value W(t0)
- `Wend`: Ending value W(tend)
- `Z0`, `Zend`: Optional auxiliary process values

# Examples
```julia
# Stock price bridge from \$100 to \$110 over 1 year
bridge = GeometricBrownianBridge(0.05, 0.2, 0.0, 1.0, 100.0, 110.0)
```
"""
function GeometricBrownianBridge(
        μ, σ, t0, tend, W0, Wend, Z0 = nothing, Zend = nothing;
        kwargs...
    )
    W = GeometricBrownianMotionProcess(μ, σ, t0, W0, Z0; kwargs...)
    h = tend - t0
    Wh = Wend - W0
    if Z0 !== nothing
        Zh = Zend - Z0
    else
        Zh = nothing
    end
    push!(W.S₁, (h, Wh, Zh))
    push!(W.reinitS₁, (h, Wh, Zh))
    return W
end

"""
    GeometricBrownianBridge!(μ, σ, t0, tend, W0, Wend,
        Z0 = nothing, Zend = nothing; kwargs...)

In-place version of [`GeometricBrownianBridge`](@ref). It conditions an
in-place geometric Brownian process on the supplied endpoints.

# Arguments
- `μ`, `σ`: Drift and diffusion coefficients.
- `t0`, `tend`: Start and end times.
- `W0`, `Wend`: Array-valued endpoint values.
- `Z0`, `Zend`: Optional auxiliary endpoint values.

# Keywords
Additional keywords are forwarded to `GeometricBrownianMotionProcess!`.

# Returns
A `NoiseProcess` with `isinplace(W) == true`.
"""
function GeometricBrownianBridge!(
        μ, σ, t0, tend, W0, Wh, Z0 = nothing, Zh = nothing;
        kwargs...
    )
    W = GeometricBrownianMotionProcess!(μ, σ, t0, W0, Z0; kwargs...)
    h = tend - t0
    Wh .-= W0
    if Z0 !== nothing
        Zh .-= Z0
    else
        Zh = nothing
    end
    push!(W.S₁, (h, Wh, Zh))
    push!(W.reinitS₁, (h, Wh, Zh))
    return W
end

"""
    CompoundPoissonBridge(rate, t0, tend, W0, Wend;
        rswm = RSWM(adaptivealg = :RSwM0), kwargs...)

A `CompoundPoissonBridge` is a compound Poisson process with pre-defined start
and end values.

This creates a jump process that is conditioned to have specific values at the beginning and end of the time interval.
The jumps are distributed to satisfy the endpoint constraint.

# Arguments
- `rate`: Jump rate (λ parameter)
- `t0`: Starting time
- `tend`: Ending time
- `W0`: Starting value W(t0)
- `Wend`: Ending value W(tend)

# Examples
```julia
# Jump process from 0 to 5 over unit time with rate 2.0
bridge = CompoundPoissonBridge(2.0, 0.0, 1.0, 0.0, 5.0)
```
"""
function CompoundPoissonBridge(rate, t0, tend, W0, Wend; rswm = RSWM(), kwargs...)
    W = CompoundPoissonProcess(rate, t0, W0; rswm = rswm, kwargs...)
    h = tend - t0
    Wh = Wend - W0
    push!(W.S₁, (h, Wh, nothing))
    push!(W.reinitS₁, (h, Wh, nothing))
    return W
end

"""
    CompoundPoissonBridge!(rate, t0, tend, W0, Wend;
        rswm = RSWM(adaptivealg = :RSwM0), kwargs...)

In-place version of [`CompoundPoissonBridge`](@ref).

# Arguments
- `rate`: Constant rate or `rate(u, p, t)` function.
- `t0`, `tend`: Start and end times.
- `W0`, `Wend`: Endpoint values; array endpoints are updated in-place.

# Keywords
- `rswm`: RSwM configuration. The default `:RSwM0` is appropriate for
  state-dependent rates.
- Additional keywords are forwarded to `CompoundPoissonProcess!`.

# Returns
A `NoiseProcess` with `isinplace(W) == true`.
"""
function CompoundPoissonBridge!(rate, t0, tend, W0, Wh; rswm = RSWM(), kwargs...)
    W = CompoundPoissonProcess!(rate, t0, W0; rswm = rswm, kwargs...)
    h = tend - t0
    Wh .-= W0
    push!(W.S₁, (h, Wh, nothing))
    push!(W.reinitS₁, (h, Wh, nothing))
    return W
end

"""
    OrnsteinUhlenbeckBridge(Θ, μ, σ, t0, tend, W0, Wend,
        Z0 = nothing; kwargs...)

An `OrnsteinUhlenbeckBridge` is an Ornstein-Uhlenbeck process with pre-defined
start and end values.

This creates a mean-reverting process that is conditioned to pass through specific values at the beginning
and end of the time interval.

# Arguments
- `Θ`: Mean reversion rate
- `μ`: Long-term mean
- `σ`: Volatility parameter
- `t0`: Starting time
- `tend`: Ending time
- `W0`: Starting value W(t0)
- `Wend`: Ending value W(tend)
- `Z0`: Optional auxiliary process value

# Examples
```julia
# Mean-reverting process from 1.0 to 0.5 over unit time
bridge = OrnsteinUhlenbeckBridge(2.0, 0.0, 0.3, 0.0, 1.0, 1.0, 0.5)
```
"""
function OrnsteinUhlenbeckBridge(Θ, μ, σ, t0, tend, W0, Wend, Z0 = nothing; kwargs...)
    ou = OrnsteinUhlenbeckProcess(Θ, μ, σ, t0, W0, Z0; kwargs...)
    h = tend - t0
    Wh = Wend .- W0
    push!(ou.S₁, (h, Wh, nothing))
    push!(ou.reinitS₁, (h, Wh, nothing))
    return ou
end

"""
    OrnsteinUhlenbeckBridge!(Θ, μ, σ, t0, tend, W0, Wend,
        Z0 = nothing; kwargs...)

In-place version of [`OrnsteinUhlenbeckBridge`](@ref). It conditions an
in-place exact OU process on the supplied endpoints.

# Arguments
- `Θ`, `μ`, `σ`: Mean-reversion rate, long-term mean, and diffusion coefficient.
- `t0`, `tend`: Start and end times.
- `W0`, `Wend`: Array-valued endpoint values.
- `Z0`: Optional auxiliary endpoint value.

# Keywords
Additional keywords are forwarded to `OrnsteinUhlenbeckProcess!`.

# Returns
A `NoiseProcess` with `isinplace(W) == true`.
"""
function OrnsteinUhlenbeckBridge!(Θ, μ, σ, t0, tend, W0, Wend, Z0 = nothing; kwargs...)
    ou = OrnsteinUhlenbeckProcess!(Θ, μ, σ, t0, W0, Z0; kwargs...)
    h = tend - t0
    Wh = Wend .- W0
    push!(ou.S₁, (h, Wh, nothing))
    push!(ou.reinitS₁, (h, Wh, nothing))
    return ou
end
