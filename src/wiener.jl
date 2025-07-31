const one_over_sqrt2 = 1 / sqrt(2)
"""
    wiener_randn(rng::AbstractRNG, ::Type{T}) where {T}

Generate a random number from the standard normal distribution for type T.

# Arguments
- `rng`: Random number generator
- `T`: Type of the random number to generate

# Returns
A random number of type T from the standard normal distribution
"""
@inline wiener_randn(rng::AbstractRNG, ::Type{T}) where {T} = randn(rng, T)

"""
    wiener_randn(rng::AbstractRNG, proto::AbstractArray{T}) where {T <: Number}

Generate an array of random numbers from the standard normal distribution, matching the size of the prototype array.

# Arguments
- `rng`: Random number generator
- `proto`: Prototype array whose size determines the output size

# Returns
An array of random numbers from the standard normal distribution with the same size as proto
"""
@inline function wiener_randn(rng::AbstractRNG, proto::AbstractArray{T}) where {T <: Number}
    randn(rng, T, size(proto))
end
"""
    wiener_randn(rng::AbstractRNG, proto::T) where {T <: StaticArraysCore.SArray}

Generate a static array of random numbers from the standard normal distribution.

# Arguments
- `rng`: Random number generator
- `proto`: Prototype static array

# Returns
A static array of the same type as proto filled with standard normal random numbers
"""
@inline function wiener_randn(rng::AbstractRNG,
        proto::T) where {T <: StaticArraysCore.SArray}
    randn(rng, T)
end
"""
    wiener_randn(rng::AbstractRNG, proto)

Generate random numbers from the standard normal distribution for arbitrary types.

# Arguments
- `rng`: Random number generator
- `proto`: Prototype object whose type and size determine the output

# Returns
Random values converted to the same type as proto
"""
@inline function wiener_randn(rng::AbstractRNG, proto)
    convert(typeof(proto), randn(rng, size(proto)))
end
"""
    wiener_randn!(rng::AbstractRNG, rand_vec::AbstractArray)

Fill an array with random numbers from the standard normal distribution in-place.

# Arguments
- `rng`: Random number generator
- `rand_vec`: Array to fill with random values

# Returns
The modified rand_vec filled with standard normal random numbers
"""
@inline wiener_randn!(rng::AbstractRNG, rand_vec::AbstractArray) = randn!(rng, rand_vec)
"""
    wiener_randn!(rng::AbstractRNG, rand_vec)

Fill an arbitrary container with random numbers from the standard normal distribution in-place using broadcasting.

# Arguments
- `rng`: Random number generator (not used in this fallback)
- `rand_vec`: Container to fill with random values

# Returns
The modified rand_vec filled with standard normal random numbers
"""
@inline function wiener_randn!(rng::AbstractRNG, rand_vec)
    rand_vec .= Base.Broadcast.Broadcasted(randn, ())
end

"""
    wiener_randn!(rng::AbstractRNG, rand_vec::GPUArraysCore.AbstractGPUArray)

Fill a GPU array with random numbers from the standard normal distribution in-place.

This specialized method works for GPUs because it doesn't pass the RNG to the GPU kernel,
which may not be supported on all GPU backends.

# Arguments
- `rng`: Random number generator (not passed to GPU)
- `rand_vec`: GPU array to fill with random values

# Returns
The modified rand_vec filled with standard normal random numbers
"""
@inline function wiener_randn!(rng::AbstractRNG, rand_vec::GPUArraysCore.AbstractGPUArray)
    randn!(rand_vec)
end

"""
    wiener_randn!(y::AbstractRNG, x::AbstractArray{<:Complex{T}}) where {T <: Number}

Fill an array of complex numbers with random values from the standard complex normal distribution in-place.

Each complex number is generated as (a + bi)/√2 where a and b are independent standard normal random variables.

# Arguments
- `y`: Random number generator
- `x`: Array of complex numbers to fill

# Returns
The modified array filled with complex normal random numbers
"""
@inline function wiener_randn!(y::AbstractRNG,
        x::AbstractArray{<:Complex{T}}) where {T <: Number}
    # Remove loop
    @inbounds for i in eachindex(x)
        x[i] = convert(T, one_over_sqrt2) * (randn(y, T) + im * randn(y, T))
    end
end

"""
    WHITE_NOISE_DIST(dW, W, dt, u, p, t, rng)

Generate white noise distributed according to N(0, dt) for use in Wiener processes.

# Arguments
- `dW`: Noise increment container
- `W`: Current noise value (unused for white noise)
- `dt`: Time step
- `u`: Current state (for state-dependent noise)
- `p`: Parameters
- `t`: Current time
- `rng`: Random number generator

# Returns
Random values distributed as N(0, dt), scaled by √dt
"""
@inline function WHITE_NOISE_DIST(dW, W, dt, u, p, t, rng)
    if dW isa AbstractArray && !(dW isa StaticArraysCore.SArray)
        return @fastmath sqrt(abs(dt)) * wiener_randn(rng, dW)
    else
        return @fastmath sqrt(abs(dt)) * wiener_randn(rng, typeof(dW))
    end
end

"""
    WHITE_NOISE_BRIDGE(dW, W, W0, Wh, q, h, u, p, t, rng)

Generate white noise for Brownian bridge interpolation between two points.

# Arguments
- `dW`: Noise increment container
- `W`: Current noise value
- `W0`: Starting noise value
- `Wh`: Target noise value at end of interval
- `q`: Interpolation parameter (0 to 1)
- `h`: Total time interval
- `u`, `p`, `t`: State, parameters, and time (for compatibility)
- `rng`: Random number generator

# Returns
Interpolated noise value that maintains correct distribution
"""
function WHITE_NOISE_BRIDGE(dW, W, W0, Wh, q, h, u, p, t, rng)
    if dW isa AbstractArray
        return @fastmath sqrt((1 - q) * q * abs(h)) * wiener_randn(rng, dW) + q * Wh
    else
        return @fastmath sqrt((1 - q) * q * abs(h)) * wiener_randn(rng, typeof(dW)) + q * Wh
    end
end

"""
    VBT_BRIDGE(dW, W, W0, Wh, q, h, u, p, t, rng)

Generate noise for Virtual Brownian Tree (VBT) bridge interpolation.

The VBT bridge is a memory-efficient method for generating Brownian paths that
can be evaluated at arbitrary time points without storing the entire path.

# Arguments
- `dW`: Noise increment container
- `W`: Current noise value
- `W0`: Starting noise value
- `Wh`: Target noise value at end of interval
- `q`: Interpolation parameter (0 to 1)
- `h`: Total time interval
- `u`, `p`, `t`: State, parameters, and time (for compatibility)
- `rng`: Random number generator

# Returns
Interpolated noise value using the VBT bridge formula
"""
function VBT_BRIDGE(dW, W, W0, Wh, q, h, u, p, t, rng)
    if dW isa AbstractArray
        return @fastmath sqrt((1 - q) * q * abs(h)) * wiener_randn(rng, dW) + q * (Wh + W0)
    else
        return @fastmath sqrt((1 - q) * q * abs(h)) * wiener_randn(rng, typeof(dW)) +
                         q * (Wh + W0)
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
function WienerProcess(t0, W0, Z0 = nothing; kwargs...)
    NoiseProcess{false}(t0, W0, Z0, WHITE_NOISE_DIST, WHITE_NOISE_BRIDGE; kwargs...)
end

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
function SimpleWienerProcess(t0, W0, Z0 = nothing; kwargs...)
    SimpleNoiseProcess{false}(t0, W0, Z0, WHITE_NOISE_DIST, WHITE_NOISE_BRIDGE; kwargs...)
end

"""
    INPLACE_WHITE_NOISE_DIST(rand_vec, W, dt, u, p, t, rng)

Generate white noise distributed according to N(0, dt) in-place.

This is the in-place version of WHITE_NOISE_DIST, modifying the provided array
rather than allocating a new one.

# Arguments
- `rand_vec`: Array to fill with noise values
- `W`: Current noise value (unused for white noise)
- `dt`: Time step
- `u`: Current state (for state-dependent noise)
- `p`: Parameters
- `t`: Current time
- `rng`: Random number generator

# Effects
Modifies rand_vec to contain random values distributed as N(0, dt)
"""
function INPLACE_WHITE_NOISE_DIST(rand_vec, W, dt, u, p, t, rng)
    wiener_randn!(rng, rand_vec)
    sqrtabsdt = @fastmath sqrt(abs(dt))
    @.. rand_vec *= sqrtabsdt
end
"""
    INPLACE_WHITE_NOISE_BRIDGE(rand_vec, W, W0, Wh, q, h, u, p, t, rng)

Generate white noise for Brownian bridge interpolation in-place.

This is the in-place version of WHITE_NOISE_BRIDGE, modifying the provided array
rather than allocating a new one.

# Arguments
- `rand_vec`: Array to fill with interpolated noise values
- `W`: Current noise value
- `W0`: Starting noise value
- `Wh`: Target noise value at end of interval
- `q`: Interpolation parameter (0 to 1)
- `h`: Total time interval
- `u`, `p`, `t`: State, parameters, and time (for compatibility)
- `rng`: Random number generator

# Effects
Modifies rand_vec to contain interpolated noise values
"""
function INPLACE_WHITE_NOISE_BRIDGE(rand_vec, W, W0, Wh, q, h, u, p, t, rng)
    wiener_randn!(rng, rand_vec)
    #rand_vec .= sqrt((1.-q).*q.*abs(h)).*rand_vec.+q.*Wh
    sqrtcoeff = @fastmath sqrt((1 - q) * q * abs(h))
    @.. rand_vec = sqrtcoeff * rand_vec + q * Wh
end

"""
    INPLACE_VBT_BRIDGE(rand_vec, W, W0, Wh, q, h, u, p, t, rng)

Generate noise for Virtual Brownian Tree (VBT) bridge interpolation in-place.

This is the in-place version of VBT_BRIDGE, modifying the provided array
rather than allocating a new one.

# Arguments
- `rand_vec`: Array to fill with interpolated noise values
- `W`: Current noise value
- `W0`: Starting noise value
- `Wh`: Target noise value at end of interval
- `q`: Interpolation parameter (0 to 1)
- `h`: Total time interval
- `u`, `p`, `t`: State, parameters, and time (for compatibility)
- `rng`: Random number generator

# Effects
Modifies rand_vec to contain VBT-interpolated noise values
"""
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
function WienerProcess!(t0, W0, Z0 = nothing; kwargs...)
    NoiseProcess{true}(t0, W0, Z0, INPLACE_WHITE_NOISE_DIST, INPLACE_WHITE_NOISE_BRIDGE;
        kwargs...)
end

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
function SimpleWienerProcess!(t0, W0, Z0 = nothing; kwargs...)
    SimpleNoiseProcess{true}(t0, W0, Z0, INPLACE_WHITE_NOISE_DIST,
        INPLACE_WHITE_NOISE_BRIDGE; kwargs...)
end

#### Real Valued Wiener Process. Ignores complex and the like
"""
    REAL_WHITE_NOISE_DIST(dW, W, dt, u, p, t, rng)

Generate real-valued white noise distributed according to N(0, dt).

Unlike WHITE_NOISE_DIST, this function always generates real-valued noise,
even if the input type would normally support complex values.

# Arguments
- `dW`: Noise increment container
- `W`: Current noise value (unused for white noise)
- `dt`: Time step
- `u`: Current state (for state-dependent noise)
- `p`: Parameters
- `t`: Current time
- `rng`: Random number generator

# Returns
Real-valued random numbers distributed as N(0, dt), scaled by √dt
"""
function REAL_WHITE_NOISE_DIST(dW, W, dt, u, p, t, rng)
    if dW isa AbstractArray
        return @fastmath sqrt(abs(dt)) * randn(rng, size(dW))
    else
        return @fastmath sqrt(abs(dt)) * randn(rng)
    end
end
"""
    REAL_WHITE_NOISE_BRIDGE(dW, W, W0, Wh, q, h, u, p, t, rng)

Generate real-valued white noise for Brownian bridge interpolation.

This function ensures the generated noise is always real-valued,
even for complex-valued endpoints.

# Arguments
- `dW`: Noise increment container
- `W`: Current noise value
- `W0`: Starting noise value
- `Wh`: Target noise value at end of interval
- `q`: Interpolation parameter (0 to 1)
- `h`: Total time interval
- `u`, `p`, `t`: State, parameters, and time (for compatibility)
- `rng`: Random number generator

# Returns
Real-valued interpolated noise value
"""
function REAL_WHITE_NOISE_BRIDGE(dW, W, W0, Wh, q, h, u, p, t, rng)
    if dW isa AbstractArray
        return @fastmath sqrt((1 - q) * q * abs(h)) * randn(rng, size(dW)) + q * Wh
    else
        return @fastmath sqrt((1 - q) * q * abs(h)) * randn(rng) + q * Wh
    end
end

@doc doc"""
The `RealWienerProcess` is a Brownian motion that is forced to be
real-valued. While the normal `WienerProcess` becomes complex valued
if `W0` is complex, this version is real valued for when you want to,
for example, solve an SDE defined by complex numbers where the noise
is in the reals.

```julia
RealWienerProcess(t0,W0,Z0=nothing;kwargs...)
RealWienerProcess!(t0,W0,Z0=nothing;kwargs...)
```
"""
function RealWienerProcess(t0, W0, Z0 = nothing; kwargs...)
    NoiseProcess{false}(t0, W0, Z0, REAL_WHITE_NOISE_DIST, REAL_WHITE_NOISE_BRIDGE;
        kwargs...)
end

"""
    REAL_INPLACE_WHITE_NOISE_DIST(rand_vec, W, dt, u, p, t, rng)

Generate real-valued white noise distributed according to N(0, dt) in-place.

This is the in-place version of REAL_WHITE_NOISE_DIST.

# Arguments
- `rand_vec`: Array to fill with noise values
- `W`: Current noise value (unused for white noise)
- `dt`: Time step
- `u`: Current state (for state-dependent noise)
- `p`: Parameters
- `t`: Current time
- `rng`: Random number generator

# Effects
Modifies rand_vec to contain real-valued random numbers distributed as N(0, dt)
"""
function REAL_INPLACE_WHITE_NOISE_DIST(rand_vec, W, dt, u, p, t, rng)
    sqabsdt = @fastmath sqrt(abs(dt))
    wiener_randn!(rng, rand_vec)
    @.. rand_vec *= sqabsdt
end
"""
    REAL_INPLACE_WHITE_NOISE_BRIDGE(rand_vec, W, W0, Wh, q, h, u, p, t, rng)

Generate real-valued white noise for Brownian bridge interpolation in-place.

This is the in-place version of REAL_WHITE_NOISE_BRIDGE.

# Arguments
- `rand_vec`: Array to fill with interpolated noise values
- `W`: Current noise value
- `W0`: Starting noise value
- `Wh`: Target noise value at end of interval
- `q`: Interpolation parameter (0 to 1)
- `h`: Total time interval
- `u`, `p`, `t`: State, parameters, and time (for compatibility)
- `rng`: Random number generator

# Effects
Modifies rand_vec to contain real-valued interpolated noise values
"""
function REAL_INPLACE_WHITE_NOISE_BRIDGE(rand_vec, W, W0, Wh, q, h, u, p, t, rng)
    wiener_randn!(rng, rand_vec)
    @.. rand_vec = @fastmath sqrt((1 - q) * q * abs(h)) * rand_vec + q * Wh
end

@doc doc"""
The `RealWienerProcess` is a Brownian motion that is forced to be
real-valued. While the normal `WienerProcess` becomes complex valued
if `W0` is complex, this version is real valued for when you want to,
for example, solve an SDE defined by complex numbers where the noise
is in the reals.

```julia
RealWienerProcess(t0,W0,Z0=nothing;kwargs...)
RealWienerProcess!(t0,W0,Z0=nothing;kwargs...)
```
"""
function RealWienerProcess!(t0, W0, Z0 = nothing; kwargs...)
    NoiseProcess{true}(t0, W0, Z0, REAL_INPLACE_WHITE_NOISE_DIST,
        REAL_INPLACE_WHITE_NOISE_BRIDGE; kwargs...)
end
