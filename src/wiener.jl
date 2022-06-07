const one_over_sqrt2 = 1 / sqrt(2)
@inline wiener_randn(rng::AbstractRNG, ::Type{T}) where {T} = randn(rng, T)
@inline function wiener_randn(rng::AbstractRNG, proto::Array{T}) where {T}
  randn(rng, size(proto))
end
@inline function wiener_randn(rng::AbstractRNG, proto::T) where {T<:SArray}
  randn(rng, T)
end
@inline function wiener_randn(rng::AbstractRNG, proto)
  convert(typeof(proto), randn(rng, size(proto)))
end
@inline wiener_randn!(rng::AbstractRNG, rand_vec::Array) = randn!(rng, rand_vec)
@inline wiener_randn(y::AbstractRNG, ::Type{Complex{T}}) where {T} = convert(T, one_over_sqrt2) * (randn(y, T) + im * randn(y, T))


@inline wiener_randn!(rng::AbstractRNG, rand_vec) = rand_vec .= Base.Broadcast.Broadcasted(randn, ())

# This fallback works for GPUs because it doesn't assume we can pass an RNG
@inline wiener_randn!(rng::AbstractRNG, rand_vec::GPUArrays.AbstractGPUArray) = randn!(rand_vec)

@inline function wiener_randn!(y::AbstractRNG, x::AbstractArray{<:Complex{T}}) where {T<:Number}
  # Remove loop
  @inbounds for i in eachindex(x)
    x[i] = convert(T, one_over_sqrt2) * (randn(y, T) + im * randn(y, T))
  end
end

@inline function WHITE_NOISE_DIST(dW, W, dt, u, p, t, rng)
  if typeof(dW) <: AbstractArray && !(typeof(dW) <: SArray)
    return @fastmath sqrt(abs(dt)) * wiener_randn(rng, dW)
  else
    return @fastmath sqrt(abs(dt)) * wiener_randn(rng, typeof(dW))
  end
end

function WHITE_NOISE_BRIDGE(dW, W, W0, Wh, q, h, u, p, t, rng)
  if typeof(dW) <: AbstractArray
    return @fastmath sqrt((1 - q) * q * abs(h)) * wiener_randn(rng, dW) + q * Wh
  else
    return @fastmath sqrt((1 - q) * q * abs(h)) * wiener_randn(rng, typeof(dW)) + q * Wh
  end
end

function VBT_BRIDGE(dW, W, W0, Wh, q, h, u, p, t, rng)
  if typeof(dW) <: AbstractArray
    return @fastmath sqrt((1 - q) * q * abs(h)) * wiener_randn(rng, dW) + q * (Wh + W0)
  else
    return @fastmath sqrt((1 - q) * q * abs(h)) * wiener_randn(rng, typeof(dW)) + q * (Wh + W0)
  end
end

@doc doc"""
The `WienerProcess`, also known as Brownian motion, or
the noise in the Langevin equation, is the stationary process with
white noise increments and a distribution `N(0,dt)`. The constructor is:

```julia
WienerProcess(t0,W0,Z0=nothing;kwargs...)
WienerProcess!(t0,W0,Z0=nothing;kwargs...)
```
"""
WienerProcess(t0, W0, Z0=nothing; kwargs...) = NoiseProcess{false}(t0, W0, Z0, WHITE_NOISE_DIST, WHITE_NOISE_BRIDGE; kwargs...)

@doc doc"""
The `SimpleWienerProcess`, also known as Brownian motion, or
the noise in the Langevin equation, is the stationary process with
white noise increments and a distribution `N(0,dt)`. The constructor is:

```julia
SimpleWienerProcess(t0,W0,Z0=nothing;kwargs...)
SimpleWienerProcess(t0,W0,Z0=nothing;kwargs...)
```

Unlike WienerProcess, this uses the SimpleNoiseProcess and thus does not
support adaptivity, but is slightly more lightweight.
"""
SimpleWienerProcess(t0, W0, Z0=nothing; kwargs...) = SimpleNoiseProcess{false}(t0, W0, Z0, WHITE_NOISE_DIST, WHITE_NOISE_BRIDGE; kwargs...)

function INPLACE_WHITE_NOISE_DIST(rand_vec, W, dt, u, p, t, rng)
  wiener_randn!(rng, rand_vec)
  sqrtabsdt = @fastmath sqrt(abs(dt))
  @.. rand_vec *= sqrtabsdt
end
function INPLACE_WHITE_NOISE_BRIDGE(rand_vec, W, W0, Wh, q, h, u, p, t, rng)
  wiener_randn!(rng, rand_vec)
  #rand_vec .= sqrt((1.-q).*q.*abs(h)).*rand_vec.+q.*Wh
  sqrtcoeff = @fastmath sqrt((1 - q) * q * abs(h))
  @.. rand_vec = sqrtcoeff * rand_vec + q * Wh
end

function INPLACE_VBT_BRIDGE(rand_vec, W, W0, Wh, q, h, u, p, t, rng)
  wiener_randn!(rng, rand_vec)
  #rand_vec .= sqrt((1.-q).*q.*abs(h)).*rand_vec.+q.*Wh
  sqrtcoeff = @fastmath sqrt((1 - q) * q * abs(h))
  @.. rand_vec = sqrtcoeff * rand_vec + q * (W0 + Wh)
end

@doc doc"""
The `WienerProcess`, also known as Brownian motion, or
the noise in the Langevin equation, is the stationary process with
white noise increments and a distribution `N(0,dt)`. The constructor is:

```julia
WienerProcess(t0,W0,Z0=nothing;kwargs...)
WienerProcess!(t0,W0,Z0=nothing;kwargs...)
```
"""
WienerProcess!(t0, W0, Z0=nothing; kwargs...) = NoiseProcess{true}(t0, W0, Z0, INPLACE_WHITE_NOISE_DIST, INPLACE_WHITE_NOISE_BRIDGE; kwargs...)

@doc doc"""
The `SimpleWienerProcess`, also known as Brownian motion, or
the noise in the Langevin equation, is the stationary process with
white noise increments and a distribution `N(0,dt)`. The constructor is:

```julia
SimpleWienerProcess(t0,W0,Z0=nothing;kwargs...)
SimpleWienerProcess(t0,W0,Z0=nothing;kwargs...)
```

Unlike WienerProcess, this uses the SimpleNoiseProcess and thus does not
support adaptivity, but is slightly more lightweight.
"""
SimpleWienerProcess!(t0, W0, Z0=nothing; kwargs...) = SimpleNoiseProcess{true}(t0, W0, Z0, INPLACE_WHITE_NOISE_DIST, INPLACE_WHITE_NOISE_BRIDGE; kwargs...)


#### Real Valued Wiener Process. Ignores complex and the like
function REAL_WHITE_NOISE_DIST(dW, W, dt, u, p, t, rng)
  if typeof(dW) <: AbstractArray
    return @fastmath sqrt(abs(dt)) * randn(rng, size(dW))
  else
    return @fastmath sqrt(abs(dt)) * randn(rng)
  end
end
function REAL_WHITE_NOISE_BRIDGE(dW, W, W0, Wh, q, h, u, p, t, rng)
  if typeof(dW) <: AbstractArray
    return @fastmath sqrt((1 - q) * q * abs(h)) * randn(rng, size(dW)) + q * Wh
  else
    return @fastmath sqrt((1 - q) * q * abs(h)) * randn(rng) + q * Wh
  end
end

@doc doc"""
The `RealWienerProcess` is a Brownian motion that is forced to be
real-valued. While the normal `WienerProcess` becomes complex valued
if `W0` is complex, this verion is real valued for when you want to,
for example, solve an SDE defined by complex numbers where the noise
is in the reals.

```julia
RealWienerProcess(t0,W0,Z0=nothing;kwargs...)
RealWienerProcess!(t0,W0,Z0=nothing;kwargs...)
```
"""
RealWienerProcess(t0, W0, Z0=nothing; kwargs...) = NoiseProcess{false}(t0, W0, Z0, REAL_WHITE_NOISE_DIST, REAL_WHITE_NOISE_BRIDGE; kwargs...)

function REAL_INPLACE_WHITE_NOISE_DIST(rand_vec, W, dt, u, p, t, rng)
  sqabsdt = @fastmath sqrt(abs(dt))
  wiener_randn!(rng, rand_vec)
  @.. rand_vec *= sqabsdt
end
function REAL_INPLACE_WHITE_NOISE_BRIDGE(rand_vec, W, W0, Wh, q, h, u, p, t, rng)
  wiener_randn!(rng, rand_vec)
  @.. rand_vec = @fastmath sqrt((1 - q) * q * abs(h)) * rand_vec + q * Wh
end

@doc doc"""
The `RealWienerProcess` is a Brownian motion that is forced to be
real-valued. While the normal `WienerProcess` becomes complex valued
if `W0` is complex, this verion is real valued for when you want to,
for example, solve an SDE defined by complex numbers where the noise
is in the reals.

```julia
RealWienerProcess(t0,W0,Z0=nothing;kwargs...)
RealWienerProcess!(t0,W0,Z0=nothing;kwargs...)
```
"""
RealWienerProcess!(t0, W0, Z0=nothing; kwargs...) = NoiseProcess{true}(t0, W0, Z0, REAL_INPLACE_WHITE_NOISE_DIST, REAL_INPLACE_WHITE_NOISE_BRIDGE; kwargs...)
