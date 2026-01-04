@doc doc"""
A `BrownianBridge` process is a Wiener process with a pre-defined start and end
value. This process is distribution exact and back be back interpolated exactly
as well. The constructor is:

```julia
BrownianBridge(t0,tend,W0,Wend,Z0=nothing,Zend=nothing;kwargs...)
BrownianBridge!(t0,tend,W0,Wend,Z0=nothing,Zend=nothing;kwargs...)
```

where `W(t0)=W₀`, `W(tend)=Wend`, and likewise for the `Z` process if defined.
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

@doc doc"""
A `BrownianBridge` process is a Wiener process with a pre-defined start and end
value. This process is distribution exact and back be back interpolated exactly
as well. The constructor is:

```julia
BrownianBridge(t0,tend,W0,Wend,Z0=nothing,Zend=nothing;kwargs...)
BrownianBridge!(t0,tend,W0,Wend,Z0=nothing,Zend=nothing;kwargs...)
```

where `W(t0)=W₀`, `W(tend)=Wend`, and likewise for the `Z` process if defined.
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

@doc doc"""
A `GeometricBrownianBridge` is a geometric Brownian motion process with pre-defined start and end values.

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
# Stock price bridge from $100 to $110 over 1 year
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

@doc doc"""
In-place version of `GeometricBrownianBridge`.

See `GeometricBrownianBridge` for details.
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

@doc doc"""
A `CompoundPoissonBridge` is a compound Poisson process with pre-defined start and end values.

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
function CompoundPoissonBridge(rate, t0, tend, W0, Wend; kwargs...)
    W = CompoundPoissonProcess(rate, t0, W0; kwargs...)
    h = tend - t0
    Wh = Wend - W0
    push!(W.S₁, (h, Wh, nothing))
    push!(W.reinitS₁, (h, Wh, nothing))
    return W
end

@doc doc"""
In-place version of `CompoundPoissonBridge`.

See `CompoundPoissonBridge` for details.
"""
function CompoundPoissonBridge!(rate, t0, tend, W0, Wh; kwargs...)
    W = CompoundPoissonProcess!(rate, t0, W0; kwargs...)
    h = tend - t0
    Wh .-= W0
    push!(W.S₁, (h, Wh, nothing))
    push!(W.reinitS₁, (h, Wh, nothing))
    return W
end

@doc doc"""
An `OrnsteinUhlenbeckBridge` is an Ornstein-Uhlenbeck process with pre-defined start and end values.

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

@doc doc"""
In-place version of `OrnsteinUhlenbeckBridge`.

See `OrnsteinUhlenbeckBridge` for details.
"""
function OrnsteinUhlenbeckBridge!(Θ, μ, σ, t0, tend, W0, Wend, Z0 = nothing; kwargs...)
    ou = OrnsteinUhlenbeckProcess!(Θ, μ, σ, t0, W0, Z0; kwargs...)
    h = tend - t0
    Wh = Wend .- W0
    push!(ou.S₁, (h, Wh, nothing))
    push!(ou.reinitS₁, (h, Wh, nothing))
    return ou
end
