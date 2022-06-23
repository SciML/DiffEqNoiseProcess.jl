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
    if Z0 != nothing
        Zh = Zend - Z0
    else
        Zh = nothing
    end
    push!(W.S₁, (h, Wh, Zh))
    W
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
    if Z0 != nothing
        Zh .-= Z0
    else
        Zh = nothing
    end
    push!(W.S₁, (h, Wh, Zh))
    W
end

function GeometricBrownianBridge(
    μ,
    σ,
    t0,
    tend,
    W0,
    Wend,
    Z0 = nothing,
    Zend = nothing;
    kwargs...,
)
    W = GeometricBrownianMotionProcess(μ, σ, t0, W0, Z0; kwargs...)
    h = tend - t0
    Wh = Wend - W0
    if Z0 != nothing
        Zh = Zend - Z0
    else
        Zh = nothing
    end
    push!(W.S₁, (h, Wh, Zh))
    W
end

function GeometricBrownianBridge!(
    μ,
    σ,
    t0,
    tend,
    W0,
    Wh,
    Z0 = nothing,
    Zh = nothing;
    kwargs...,
)
    W = GeometricBrownianMotionProcess!(μ, σ, t0, W0, Z0; kwargs...)
    h = tend - t0
    Wh .-= W0
    if Z0 != nothing
        Zh .-= Z0
    else
        Zh = nothing
    end
    push!(W.S₁, (h, Wh, Zh))
    W
end

function CompoundPoissonBridge(rate, t0, tend, W0, Wend; kwargs...)
    W = CompoundPoissonProcess(rate, t0, W0; kwargs...)
    h = tend - t0
    Wh = Wend - W0
    push!(W.S₁, (h, Wh, nothing))
    W
end

function CompoundPoissonBridge!(rate, t0, tend, W0, Wh; kwargs...)
    W = CompoundPoissonProcess!(rate, t0, W0; kwargs...)
    h = tend - t0
    Wh .-= W0
    push!(W.S₁, (h, Wh, nothing))
    W
end
