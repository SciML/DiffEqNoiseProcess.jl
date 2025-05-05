isinplace(W::AbstractNoiseProcess{T, N, S, inplace}) where {T, N, S, inplace} = inplace

"""
```julia
NoiseProcess{T, N, Tt, T2, T3, ZType, F, F2, inplace, S1, S2, RSWM, C, RNGType} <:
{T, N, Vector{T2}, inplace}
```

A `NoiseProcess` is a type defined as:

```julia
NoiseProcess(t0, W0, Z0, dist, bridge;
    iip = SciMLBase.isinplace(dist, 3),
    rswm = RSWM(), save_everystep = true,
    rng = Xorshifts.Xoroshiro128Plus(rand(UInt64)),
    reset = true, reseed = true)
```

  - `t0` is the first timepoint
  - `W0` is the first value of the process.
  - `Z0` is the first value of the pseudo-process. This is necessary for higher
    order algorithms. If it's not needed, set to `nothing`.
  - `dist` the distribution of the steps over time.
  - `bridge` the bridging distribution. Optional, but required for adaptivity and interpolating
    at new values.
  - `covariance` is the covariance matrix of the noise process. If not provided, the noise
    is assumed to be uncorrelated in each variable.
  - `save_everystep` whether to save every step of the Brownian timeseries.
  - `rng` the local RNG used for generating the random numbers.
  - `reset` whether to reset the process with each solve.
  - `reseed` whether to reseed the process with each solve.

The signature for the `dist` is:

```julia
dist!(rand_vec, dW, W, dt, u, p, t, rng)
```

for inplace functions, and:

```julia
rand_vec = dist(dW, W, dt, u, p, t, rng)
```

otherwise. The signature for `bridge` is:

```julia
bridge!(rand_vec, dW, W, W0, Wh, q, h, u, p, t, rng)
```

and the out of place syntax is:

```julia
rand_vec = bridge(dW, W, W0, Wh, q, h, u, p, t, rng)
```

Here, `W` is the noise process, `W0` is the left side of the current interval,
`Wh` is the right side of the current interval, `h` is the interval length,
and `q` is the proportion from the left where the interpolation is occurring.

## Direct Construction Example

The easiest way to show how to directly construct a `NoiseProcess` is by example.
Here we will show how to directly construct a `NoiseProcess` which generates
Gaussian white noise.

This is the noise process, that uses `randn!`. A special dispatch is added for
complex numbers for `(randn()+im*randn())/sqrt(2)`. This function is
`DiffEqNoiseProcess.wiener_randn` (or with `!` respectively).

The first function that must be defined is the noise distribution. This is how
to generate ``W(t+dt)`` given that we know ``W(x)`` for ``x∈[t₀,t]``. For Gaussian
white noise, we know that

```math
W(dt) ∼ N(0,dt)
```

for ``W(0)=0`` which defines the stepping distribution. Thus, its noise distribution
function is:

```julia
@inline function WHITE_NOISE_DIST(dW, W, dt, u, p, t, rng)
    if W.dW isa AbstractArray && !(W.dW isa SArray)
        return @fastmath sqrt(abs(dt)) * wiener_randn(rng, W.dW)
    else
        return @fastmath sqrt(abs(dt)) * wiener_randn(rng, typeof(W.dW))
    end
end
```

for the out of place versions, and for the inplace versions

```julia
function INPLACE_WHITE_NOISE_DIST(rand_vec, dW, W, dt, u, p, t, rng)
    wiener_randn!(rng, rand_vec)
    sqrtabsdt = @fastmath sqrt(abs(dt))
    @. rand_vec *= sqrtabsdt
end
```

Optionally, we can provide a bridging distribution. This is the distribution of
``W(qh)`` for ``q∈[0,1]`` given that we know ``W(0)=0`` and ``W(h)=Wₕ``. For
Brownian motion, this is known as the Brownian Bridge, and is well known to have
the distribution:

```math
W(qh) ∼ N(qWₕ,(1-q)qh)
```

Thus, we have the out-of-place and in-place versions as:

```julia
function WHITE_NOISE_BRIDGE(dW, W, W0, Wh, q, h, u, p, t, rng)
    if W.dW isa AbstractArray
        return @fastmath sqrt((1 - q) * q * abs(h)) * wiener_randn(rng, W.dW) + q * Wh
    else
        return @fastmath sqrt((1 - q) * q * abs(h)) * wiener_randn(rng, typeof(W.dW)) +
                         q * Wh
    end
end
function INPLACE_WHITE_NOISE_BRIDGE(rand_vec, dW, W, W0, Wh, q, h, u, p, t, rng)
    wiener_randn!(rng, rand_vec)
    #rand_vec .= sqrt((1.-q).*q.*abs(h)).*rand_vec.+q.*Wh
    sqrtcoeff = @fastmath sqrt((1 - q) * q * abs(h))
    @. rand_vec = sqrtcoeff * rand_vec + q * Wh
end
```

These functions are then placed in a noise process:

```julia
NoiseProcess(t0, W0, Z0, WHITE_NOISE_DIST, WHITE_NOISE_BRIDGE; kwargs)
NoiseProcess(t0, W0, Z0, INPLACE_WHITE_NOISE_DIST, INPLACE_WHITE_NOISE_BRIDGE; kwargs)
```

Notice that we can optionally provide an alternative adaptive algorithm for the
timestepping rejections. `RSWM()` defaults to the Rejection Sampling with Memory
3 algorithm (RSwM3).

Note that the standard constructors are simply:

```julia
function WienerProcess(t0, W0, Z0 = nothing)
    NoiseProcess(t0, W0, Z0, WHITE_NOISE_DIST, WHITE_NOISE_BRIDGE; kwargs)
end
function WienerProcess!(t0, W0, Z0 = nothing)
    NoiseProcess(t0, W0, Z0, INPLACE_WHITE_NOISE_DIST, INPLACE_WHITE_NOISE_BRIDGE; kwargs)
end
```

These will generate a Wiener process, which can be stepped with `step!(W,dt)`, and interpolated as `W(t)`.
"""
mutable struct NoiseProcess{
    T, N, Tt, T2, T3, ZType, F, F2, CovType, inplace, S1, S2, RSWM, C,
    RNGType} <: AbstractNoiseProcess{T, N, Vector{T2}, inplace}
    dist::F
    bridge::F2
    t::Vector{Tt}
    u::Vector{T2} # Aliased pointer to W for the AbstractVectorOfArray interface
    W::Vector{T2}
    Z::ZType
    curt::Tt
    curW::T2
    curZ::T3
    dt::Tt
    dW::T2
    dZ::T3
    dWtilde::T2
    dZtilde::T3
    dWtmp::T2
    dZtmp::T3
    covariance::CovType
    S₁::S1
    S₂::S2
    reinitS₁::S1
    rswm::RSWM
    maxstacksize::Int
    maxstacksize2::Int
    save_everystep::Bool
    iter::Int
    rng::RNGType
    reset::Bool
    reseed::Bool
    continuous::Bool
    cache::C
end

function NoiseProcess{iip}(t0, W0, Z0, dist, bridge;
        rswm = RSWM(), save_everystep = true, covariance = nothing,
        rng = Xorshifts.Xoroshiro128Plus(rand(UInt64)),
        reset = true, reseed = true, continuous = true,
        cache = nothing) where {iip}
    S₁ = ResettableStacks.ResettableStack{iip}(Tuple{typeof(t0), typeof(W0), typeof(Z0)
    })
    S₂ = ResettableStacks.ResettableStack{iip}(Tuple{typeof(t0), typeof(W0), typeof(Z0)
    })
    reinitS₁ = ResettableStacks.ResettableStack{iip}(Tuple{typeof(t0), typeof(W0),
        typeof(Z0)
    })
    if Z0 == nothing
        Z = nothing
        curZ = nothing
        dZ = nothing
        dZtilde = nothing
        dZtmp = nothing
    else
        Z = [copy(Z0)]
        curZ = copy(Z0)
        dZ = copy(Z0)
        dZtilde = copy(Z0)
        dZtmp = copy(Z0)
    end
    W = [copy(W0)]
    N = length((size(W0)..., length(W)))
    NoiseProcess{eltype(eltype(W0)), N, typeof(t0), typeof(W0), typeof(dZ), typeof(Z),
        typeof(dist), typeof(bridge), typeof(covariance),
        iip, typeof(S₁), typeof(S₂), typeof(rswm), typeof(cache), typeof(rng)}(dist,
        bridge,
        [
            t0
        ],
        W,
        W,
        Z,
        t0,
        copy(W0),
        curZ,
        t0,
        copy(W0),
        dZ,
        copy(W0),
        dZtilde,
        copy(W0),
        dZtmp,
        covariance,
        S₁,
        S₂,
        reinitS₁,
        rswm,
        0,
        0,
        save_everystep,
        0,
        rng,
        reset,
        reseed,
        continuous,
        cache)
end

function vec_NoiseProcess(W::NoiseProcess{T, N, Tt}) where {T, N, Tt}
    iip = DiffEqBase.isinplace(W)
    _W = vec.(W.u)
    Wtype = typeof(vec(W.W[1]))

    if W.curZ === nothing
        Ztype = typeof(nothing)
        Z = nothing
        curZ = nothing
        dZ = nothing
        dZtilde = nothing
        dZtmp = nothing
    else
        Ztype = typeof(vec(W.Z[1]))
        Z = vec.(W.Z)
        curZ = vec(W.curZ)
        dZ = vec(W.dZ)
        dZtilde = vec(W.dZtilde)
        dZtmp = vec(W.dZtmp)
    end

    S₁ = ResettableStacks.ResettableStack{iip}(Tuple{typeof(W.curt), Wtype, Ztype})
    S₂ = ResettableStacks.ResettableStack{iip}(Tuple{typeof(W.curt), Wtype, Ztype})
    reinitS₁ = ResettableStacks.ResettableStack{iip}(Tuple{typeof(W.curt), Wtype, Ztype})

    NoiseProcess{T, N, Tt, Wtype, typeof(dZ), typeof(Z),
        typeof(W.dist), typeof(W.bridge), typeof(W.covariance),
        iip, typeof(S₁), typeof(S₂), typeof(W.rswm), typeof(W.cache), typeof(W.rng)
    }(W.dist,
        W.bridge,
        W.t, _W,
        _W, Z, W.curt,
        vec(W.curW),
        curZ, W.curt,
        vec(W.dW),
        dZ,
        vec(W.dWtilde),
        dZtilde,
        vec(W.dWtmp),
        dZtmp,
        W.covariance,
        S₁, S₂,
        reinitS₁,
        W.rswm, W.maxstacksize,
        W.maxstacksize2,
        W.save_everystep,
        W.iter, W.rng,
        W.reset,
        W.reseed,
        W.continuous,
        W.cache)
end

(W::NoiseProcess)(t) = interpolate!(W, nothing, nothing, t)
(W::NoiseProcess)(u, p, t) = interpolate!(W, u, p, t)
(W::NoiseProcess)(out1, out2, u, p, t) = interpolate!(out1, out2, W, u, p, t)
adaptive_alg(W::NoiseProcess) = adaptive_alg(W.rswm)

function NoiseProcess(t0, W0, Z0, dist, bridge; kwargs...)
    iip = DiffEqBase.isinplace(dist, 7)
    NoiseProcess{iip}(t0, W0, Z0, dist, bridge; kwargs...)
end

"""
```julia
SimpleNoiseProcess{T, N, Tt, T2, T3, ZType, F, F2, inplace, RNGType} <:
AbstractNoiseProcess{T, N, Vector{T2}, inplace}
```

Like `NoiseProcess` but without support for adaptivity. This makes it lightweight and slightly faster.

!!! warn

    `SimpleNoiseProcess` should not be used with adaptive SDE solvers as it will lead to incorrect results.

```julia
SimpleNoiseProcess{iip}(t0, W0, Z0, dist, bridge;
    save_everystep = true,
    rng = Xorshifts.Xoroshiro128Plus(rand(UInt64)),
    reset = true, reseed = true) where {iip}
```

  - `t0` is the first timepoint
  - `W0` is the first value of the process.
  - `Z0` is the first value of the pseudo-process. This is necessary for higher
    order algorithms. If it's not needed, set to `nothing`.
  - `dist` the distribution for the steps over time.
  - `bridge` the bridging distribution. Optional, but required for adaptivity and interpolating
    at new values.
  - `save_everystep` whether to save every step of the Brownian timeseries.
  - `rng` the local RNG used for generating the random numbers.
  - `reset` whether to reset the process with each solve.
  - `reseed` whether to reseed the process with each solve.
"""
mutable struct SimpleNoiseProcess{
    T, N, Tt, T2, T3, ZType, F, F2, CovType, inplace, RNGType} <:
               AbstractNoiseProcess{T, N, Vector{T2}, inplace}
    dist::F
    bridge::F2
    t::Vector{Tt}
    u::Vector{T2} # Aliased pointer to W for the AbstractVectorOfArray interface
    W::Vector{T2}
    Z::ZType
    curt::Tt
    curW::T2
    curZ::T3
    dt::Tt
    dW::T2
    dZ::T3
    dWtilde::T2
    dZtilde::T3
    dWtmp::T2
    dZtmp::T3
    covariance::CovType
    save_everystep::Bool
    iter::Int
    rng::RNGType
    reset::Bool
    reseed::Bool

    function SimpleNoiseProcess{iip}(t0, W0, Z0, dist, bridge;
            save_everystep = true, covariance = nothing,
            rng = Xorshifts.Xoroshiro128Plus(rand(UInt64)),
            reset = true, reseed = true) where {iip}
        if Z0 == nothing
            Z = nothing
            curZ = nothing
            dZ = nothing
            dZtilde = nothing
            dZtmp = nothing
        else
            Z = [copy(Z0)]
            curZ = copy(Z0)
            dZ = copy(Z0)
            dZtilde = copy(Z0)
            dZtmp = copy(Z0)
        end
        W = [copy(W0)]
        N = length((size(W0)..., length(W)))
        new{eltype(eltype(W0)), N, typeof(t0), typeof(W0), typeof(dZ), typeof(Z),
            typeof(dist), typeof(bridge), typeof(covariance), iip, typeof(rng)}(
            dist, bridge, [t0], W, W, Z, t0,
            copy(W0), curZ, t0, copy(W0),
            dZ, copy(W0), dZtilde, copy(W0),
            dZtmp, covariance,
            save_everystep, 0, rng, reset,
            reseed)
    end
end
(W::SimpleNoiseProcess)(t) = interpolate!(W, nothing, nothing, t)
(W::SimpleNoiseProcess)(u, p, t) = interpolate!(W, u, p, t)
(W::SimpleNoiseProcess)(out1, out2, u, p, t) = interpolate!(out1, out2, W, u, p, t)

function SimpleNoiseProcess(t0, W0, Z0, dist, bridge; kwargs...)
    iip = DiffEqBase.isinplace(dist, 7)
    SimpleNoiseProcess{iip}(t0, W0, Z0, dist, bridge; kwargs...)
end

"""
```julia
NoiseWrapper{T, N, Tt, T2, T3, T4, ZType, inplace} <:
AbstractNoiseProcess{T, N, Vector{T2}, inplace}
```

This produces a new noise process from an old one, which will use its interpolation
to generate the noise. This allows you to reuse a previous noise process not just
with the same timesteps, but also with new (adaptive) timesteps as well. Thus
this is very good for doing Multi-level Monte Carlo schemes and strong
convergence testing.

## Constructor

```julia
NoiseWrapper(source::AbstractNoiseProcess{T, N, Vector{T2}, inplace};
    reset = true, reverse = false, indx = nothing) where {T, N, T2, inplace}
```

## NoiseWrapper Example

In this example, we will solve an SDE three times:

  - First, to generate a noise process
  - Second, with the same timesteps to show the values are the same
  - Third, with half-sized timesteps

First, we will generate a noise process by solving an SDE:

```julia
using StochasticDiffEq, DiffEqNoiseProcess
f1(u, p, t) = 1.01u
g1(u, p, t) = 1.01u
dt = 1 // 2^(4)
prob1 = SDEProblem(f1, g1, 1.0, (0.0, 1.0))
sol1 = solve(prob1, EM(), dt = dt, save_noise = true)
```

Now we wrap the noise into a NoiseWrapper and solve the same problem:

```julia
W2 = NoiseWrapper(sol1.W)
prob1 = SDEProblem(f1, g1, 1.0, (0.0, 1.0), noise = W2)
sol2 = solve(prob1, EM(), dt = dt)
```

We can test

```julia
@test sol1.u ≈ sol2.u
```

to see that the values are essentially equal. Now we can use the same process
to solve the same trajectory with a smaller `dt`:

```julia
W3 = NoiseWrapper(sol1.W)
prob2 = SDEProblem(f1, g1, 1.0, (0.0, 1.0), noise = W3)

dt = 1 // 2^(5)
sol3 = solve(prob2, EM(), dt = dt)
```

We can plot the results to see what this looks like:

```julia
using Plots
plot(sol1)
plot!(sol2)
plot!(sol3)
```

![noise_process](assets/noise_process.png)

In this plot, `sol2` covers up `sol1` because they hit essentially the same
values. You can see that `sol3` is similar to the others, because it's
using the same underlying noise process, just sampled much finer.

To double-check, we see that:

```julia
plot(sol1.W)
plot!(sol2.W)
plot!(sol3.W)
```

![coupled_wiener](assets/coupled_wiener.png)

the coupled Wiener processes coincide at every other time point, and the intermediate
timepoints were calculated according to a Brownian bridge.

### Adaptive NoiseWrapper Example

Here we will show that the same noise can be used with the adaptive methods
using the `NoiseWrapper`. `SRI` and `SRIW1` use slightly different error
estimators, and thus have slightly different stepping behavior. We can
see how they solve the same 2D SDE differently by using the noise
wrapper:

```julia
prob = SDEProblem(f1, g1, ones(2), (0.0, 1.0))
sol4 = solve(prob, SRI(), abstol = 1e-8, save_noise = true)

W2 = NoiseWrapper(sol4.W)
prob2 = SDEProblem(f1, g1, ones(2), (0.0, 1.0), noise = W2)
sol5 = solve(prob2, SRIW1(), abstol = 1e-8)

using Plots
plot(sol4)
plot!(sol5)
```

![SRI_SRIW1_diff](assets/SRI_SRIW1_diff.png)
"""
mutable struct NoiseWrapper{T, N, Tt, T2, T3, T4, ZType, inplace} <:
               AbstractNoiseProcess{T, N, Vector{T2}, inplace}
    t::Vector{Tt}
    u::Vector{T2}
    W::Vector{T2}
    Z::ZType
    curt::Tt
    curW::T2
    curZ::T3
    dt::Tt
    dW::T2
    dZ::T3
    source::T4
    reset::Bool
    reverse::Bool
end

function NoiseWrapper(source::AbstractNoiseProcess{T, N, Vector{T2}, inplace};
        reset = true, reverse = false,
        indx = nothing) where {T, N, T2, inplace}
    if indx === nothing
        if reverse
            indx = length(source.t)
        else
            indx = 1
        end
    end

    if source.Z == nothing
        Z = nothing
        curZ = nothing
        dZ = nothing
    else
        Z = [copy(source.Z[indx])]
        curZ = copy(source.Z[indx])
        dZ = copy(source.Z[indx])
    end
    W = [copy(source.W[indx])]

    NoiseWrapper{
        T, N, typeof(source.t[1]), typeof(source.W[1]), typeof(dZ), typeof(source),
        typeof(Z), inplace}([source.t[indx]], W, W, Z, source.t[indx],
        copy(source.W[indx]), curZ, source.t[indx],
        copy(source.W[indx]), dZ, source, reset, reverse)
end

(W::NoiseWrapper)(t) = interpolate!(W, nothing, nothing, t)
(W::NoiseWrapper)(u, p, t) = interpolate!(W, u, p, t)
function (W::NoiseWrapper)(out1, out2, u, p, t)
    interpolate!(out1, out2, W, u, p, t, reverse = W.reverse)
end
adaptive_alg(W::NoiseWrapper) = adaptive_alg(W.source)

"""
```julia
NoiseFunction{T, N, wType, zType, Tt, T2, T3, inplace} <:
AbstractNoiseProcess{T, N, nothing, inplace}
```

This allows you to use any arbitrary function `W(t)` as a `NoiseProcess`. This
will use the function lazily, only caching values required to minimize function
calls, but not storing the entire noise array. This requires an initial time point
`t0` in the domain of `W`. A second function is needed if the desired SDE algorithm
requires multiple processes.

```julia
NoiseFunction{iip}(t0, W, Z = nothing;
    noise_prototype = W(nothing, nothing, t0),
    reset = true) where {iip}
```

Additionally, one can use an in-place function `W(out1,out2,t)` for more efficient
generation of the arrays for mulitidimensional processes. When the in-place version
is used without a dispatch for the out-of-place version, the `noise_prototype`
needs to be set.

## NoiseFunction Example

The `NoiseFunction` is pretty simple: pass a function. As a silly example, we
can use `exp` as a noise process by doing:

```julia
f(u, p, t) = exp(t)
W = NoiseFunction(0.0, f)
```

If it's mulitidimensional and an in-place function is used, the `noise_prototype`
must be given. For example:

```julia
f(out, u, p, t) = (out .= exp(t))
W = NoiseFunction(0.0, f, noise_prototype = rand(4))
```

This allows you to put arbitrarily weird noise into SDEs and RODEs. Have fun.
"""
mutable struct NoiseFunction{T, N, wType, zType, Tt, T2, T3, inplace} <:
               AbstractNoiseProcess{T, N, nothing, inplace}
    W::wType
    Z::zType
    curt::Tt
    curW::T2
    curZ::T3
    dt::Tt
    dW::T2
    dZ::T3
    t0::Tt
    reset::Bool

    function NoiseFunction{iip}(t0, W, Z = nothing;
            noise_prototype = W(nothing, nothing, t0),
            reset = true) where {iip}
        curt = t0
        dt = t0
        curW = copy(noise_prototype)
        dW = copy(noise_prototype)
        if Z == nothing
            curZ = nothing
            dZ = nothing
        else
            curZ = copy(noise_prototype)
            dZ = copy(noise_prototype)
        end
        new{typeof(noise_prototype), ndims(noise_prototype), typeof(W), typeof(Z),
            typeof(curt), typeof(curW), typeof(curZ), iip}(W, Z, curt, curW, curZ,
            dt, dW, dZ, t0, reset)
    end
end

(W::NoiseFunction)(t) = W(nothing, nothing, t)
function (W::NoiseFunction)(u, p, t)
    if W.Z !== nothing
        if isinplace(W)
            out2 = similar(W.dZ)
            W.Z(out2, u, p, t)
        else
            out2 = W.Z(u, p, t)
        end
    else
        out2 = nothing
    end
    if isinplace(W)
        out1 = similar(W.dW)
        W.W(out1, u, p, t)
    else
        out1 = W.W(u, p, t)
    end
    out1, out2
end
function (W::NoiseFunction)(out1, out2, u, p, t)
    W.W(out1, u, p, t)
    W.Z !== nothing && W.Z(out2, u, p, t)
end

function NoiseFunction(t0, W, Z = nothing; kwargs...)
    iip = DiffEqBase.isinplace(W, 4)
    NoiseFunction{iip}(t0, W, Z; kwargs...)
end

"""
```julia
NoiseTransport{T, N, wType, zType, Tt, T2, T3, TRV, Trv, RNGType, inplace} <:
AbstractNoiseProcess{T, N, nothing, inplace}
```

This allows you to define stochastic processes of the form `W(t) = f(u, p, t, RV)`, where `f` is a function and `RV` represents a random variable.
This will use the function lazily, only caching values required to minimize function
calls, but not storing the entire noise array. This requires an initial time point
`t0` in the domain of `W`. A second function is needed if the desired SDE algorithm
requires multiple processes.

```julia
NoiseTransport{iip}(t0,
    W,
    RV,
    rv,
    Z = nothing;
    rng = Xorshifts.Xoroshiro128Plus(rand(UInt64)),
    reset = true,
    reseed = true,
    noise_prototype = W(nothing, nothing, t0, rv)) where {iip}
```

```julia
NoiseTransport(t0,
    W,
    RV;
    rng = Xorshifts.Xoroshiro128Plus(rand(UInt64)),
    reset = true,
    reseed = true,
    kwargs...)
```

Additionally, one can use an in-place function `W(out, u, p, t, rv)` for more efficient
generation of the arrays for mulitidimensional processes. When the in-place version
is used without a dispatch for the out-of-place version, the `noise_prototype`
needs to be set.

## NoiseTransport Example

The `NoiseTransport` requires you to pass an initial time, a transport function, and a random variable. The random variable can be either out-of-place or in-place. It is assumed it is out-of-place when the realization is a subtype of `Number`, and in-place, when it is a subtype of `AbstractArray`. Here, a random variable is any function that accepts a random number generator, in the out-of-place case (e.g. `rand(rng)`), or a random number generator and a realization to be mutated (e.g. `rand!(rng, rv)`).

An optional realization `rv` may be given. The realization `rv` is used in the first time an `AbstractRODEProblem` is solved. Subsequent runs of the same problem will draw a different realization from the random variable `RV`, unless `reseed` is set to false. In the case of a `NoiseProblem`, however, a new realization will happen at the first run already, and, in this case, `rv` can be regarded as a realization prototype, which is necessary in the case of a random vector.

As a first example, let us implement the Gaussian noise `W(t) = sin(Yt)`, where `Y` is a normal random variable.

```julia
f(u, p, t, rv) = sin(rv * t)
t0 = 0.0
W = NoiseTransport(t0, f, randn)
```

If we want to build a scalar random process out of a random vector, then an in-place version of the random vector is required, as follows. We can also use parameters in the transport function, in which case the `noise_prototype` must be given.

```julia
using Random: randn!
f(u, p, t, rv) = sin(p[1] * t + rv[1]) + cos(p[2] * t + rv[2])
t0 = 0.0
rv = randn(2)
p = (π, 2π)
W = NoiseTransport(t0, f, randn!, rv, noise_prototype = f(nothing, p, t0, rv))
```

If the random process is expected to be mulitidimensional, it is preferable to use an in-place transport function, and, in this case, the `noise_prototype` must be given. Here is an example with a scalar random vector with a beta distribution, from `Distributions.jl`.

```julia
f!(out, u, p, t, rv) = (out .= sin.(rv * t))
t0 = 0.0
RV(rng) = rand(rng, Beta(2, 3))
rv = 0.0
W = NoiseTransport(t0, f!, RV, rv, noise_prototype = zeros(4))
```

We can also have a random vector with a mulitidimensional process, in which case an in-place version of `RV` is required. For example.

```julia
using Random: randn!

function f!(out, u, p, t, v)
    out[1] = sin(v[1] * t)
    out[2] = sin(t + v[2])
    out[3] = cos(t) * v[1] + sin(t) * v[2]
    nothing
end

t0 = 0.0
RV!(rng, v) = (v[1] = randn(rng); v[2] = rand(rng))
rv = zeros(2)

W = NoiseTransport(t0, f!, RV!, rv, noise_prototype = zeros(3))
```

A `NoiseTransport` can be used as driving noise for SDEs and RODEs. Have fun!
"""
mutable struct NoiseTransport{T, N, wType, zType, Tt, T2, T3, TRV, Trv, RNGType, inplace} <:
               AbstractNoiseProcess{T, N, nothing, inplace}
    W::wType
    Z::zType
    curt::Tt
    curW::T2
    curZ::T3
    dt::Tt
    dW::T2
    dZ::T3
    t0::Tt
    RV::TRV
    rv::Trv
    rng::RNGType
    reset::Bool
    reseed::Bool

    function NoiseTransport{iip}(t0, W, RV, rv, Z = nothing;
            rng = Xorshifts.Xoroshiro128Plus(rand(UInt64)),
            reset = true, reseed = true,
            noise_prototype = W(nothing, nothing, t0, rv)) where {iip}
        curt = t0
        dt = t0
        curW = copy(noise_prototype)
        dW = copy(noise_prototype)
        if Z === nothing
            curZ = nothing
            dZ = nothing
        else
            curZ = copy(noise_prototype)
            dZ = copy(noise_prototype)
        end

        new{typeof(noise_prototype), ndims(noise_prototype), typeof(W), typeof(Z),
            typeof(curt), typeof(curW), typeof(curZ), typeof(RV), typeof(rv), typeof(rng),
            iip}(W, Z, curt, curW, curZ,
            dt, dW, dZ, t0, RV, rv, rng, reset, reseed)
    end
end

(W::NoiseTransport)(t) = W(nothing, nothing, t, W.rv)
function (W::NoiseTransport)(u, p, t, rv)
    if W.Z !== nothing
        if isinplace(W)
            out2 = similar(W.dZ)
            W.Z(out2, u, p, t, rv)
        else
            out2 = W.Z(u, p, t, rv)
        end
    else
        out2 = nothing
    end
    if isinplace(W)
        out1 = similar(W.dW)
        W.W(out1, u, p, t, rv)
    else
        out1 = W.W(u, p, t, rv)
    end
    out1, out2
end
function (W::NoiseTransport)(out1, out2, u, p, t, rv)
    W.W(out1, u, p, t, rv)
    W.Z !== nothing && W.Z(out2, u, p, t, rv)
end

function NoiseTransport(t0, W, RV, rv, Z = nothing;
        rng = Xorshifts.Xoroshiro128Plus(rand(UInt64)), reset = true,
        reseed = true, kwargs...)
    iip = DiffEqBase.isinplace(W, 5)
    NoiseTransport{iip}(t0, W, RV, rv, Z; rng, reset, reseed, kwargs...)
end

function NoiseTransport(t0, W, RV; rng = Xorshifts.Xoroshiro128Plus(rand(UInt64)),
        reset = true, reseed = true, kwargs...)
    iip = DiffEqBase.isinplace(W, 5)
    rv = RV(rng)
    Z = nothing
    NoiseTransport{iip}(t0, W, RV, rv, Z; rng, reset, reseed, kwargs...)
end

"""
A noise grid builds a noise process from arrays of points. For example, you
can generate your desired noise process as an array `W` with timepoints `t`,
and use the constructor:

```julia
NoiseGrid(t, W, Z = nothing; reset = true)
```

to build the associated noise process. This process comes with a linear
interpolation of the given points, and thus the grid does not have to match
the grid of integration. Thus, this can be used for adaptive solutions as
well. However, one must take note that the fidelity of the noise process
is linked to how fine the noise grid is determined: if the noise grid is
sparse on points compared to the integration, then its distributional
properties may be slightly perturbed by the linear interpolation. Thus, it's
suggested that the grid size at least approximates the number of
time steps in the integration to ensure accuracy.

For a one-dimensional process, `W` should be an `AbstractVector` of `Number`s.
For mulitidimensional processes, `W` should be an `AbstractVector` of the
`noise_prototype`.

## NoiseGrid

In this example, we will show you how to define your own version of Brownian
motion using an array of pre-calculated points. In normal usage, you should use
`WienerProcess` instead, since this will have distributionally-exact interpolations
while the noise grid uses linear interpolations, but this is a nice example
of the workflow.

To define a `NoiseGrid` you need to have a set of time points and a set of
values for the process. Let's define a Brownian motion on `(0.0,1.0)` with
a `dt=0.001`. To do this,

```julia
dt = 0.001
t = 0:dt:1
brownian_values = cumsum([0; [sqrt(dt) * randn() for i in 1:(length(t) - 1)]])
```

Now we build the `NoiseGrid` using these values:

```julia
W = NoiseGrid(t, brownian_values)
```

We can then pass `W` as the `noise` argument of an `SDEProblem` to use it in
an SDE.
"""
mutable struct NoiseGrid{T, N, Tt, T2, T3, ZType, RefType, inplace} <:
               AbstractNoiseProcess{T, N, Vector{T2}, inplace}
    t::Vector{Tt}
    u::Vector{T2}
    W::Vector{T2}
    Z::ZType
    curt::Tt
    curW::T2
    curZ::T3
    dt::Tt
    dW::T2
    dZ::T3
    step_setup::Bool
    reset::Bool
    cur_time::RefType
end

function NoiseGrid(t, W, Z = nothing; reset = true)
    val = W[1]
    curt = t[1]
    dt = t[1]
    curW = copy(val)
    dW = copy(val)
    if Z == nothing
        curZ = nothing
        dZ = nothing
    else
        curZ = copy(Z[1])
        dZ = copy(Z[1])
    end

    if sign(t[end] - t[1]) == 1
        cur_time = Ref(1)
    else
        cur_time = Ref(length(t))
    end

    (val isa AbstractArray && !(val isa SArray)) ? iip = true : iip = false
    NoiseGrid{typeof(val), ndims(val), typeof(dt), typeof(dW), typeof(dZ), typeof(Z),
        typeof(cur_time), iip}(t,
        W,
        W,
        Z,
        curt,
        curW,
        curZ,
        dt,
        dW,
        dZ,
        true,
        reset,
        cur_time)
end

(W::NoiseGrid)(t) = interpolate!(W, t)
(W::NoiseGrid)(u, p, t) = interpolate!(W, t)
(W::NoiseGrid)(out1, out2, u, p, t) = interpolate!(out1, out2, W, t)

"""
In many cases, one would like to define a noise process directly by a stochastic
differential equation which does not have an analytical solution. Of course,
this will not be distributionally-exact and how well the properties match
depends on how well the differential equation is integrated, but in many cases
, this can be used as a good approximation when other methods are much
more difficult.

A `NoiseApproximation` is defined by a `DEIntegrator`. The constructor for a
`NoiseApproximation` is:

```julia
NoiseApproximation(source1::DEIntegrator,
    source2::Union{DEIntegrator, Nothing} = nothing;
    reset = true)
```

The `DEIntegrator` should have a final time point of integration far enough away, such
that it will not halt during the integration. For ease of use, you can use a
final time point as `Inf`. Note that the time points do not have to match the
time points of the future integration, since the interpolant of the SDE solution
will be used. Thus, the limiting factor is error tolerance, and not hitting specific
points.

## NoiseApproximation Example

In this example, we will show how to use the `NoiseApproximation` to
build our own Geometric Brownian Motion from its stochastic differential
equation definition. In normal usage, you should use the `GeometricBrownianMotionProcess`
instead since that is more efficient and distributionally-exact.

First, let's define the `SDEProblem`. Here, we use a timespan `(0.0,Inf)` so
that the noise can be used over an indefinite integral.

```julia
const μ = 1.5
const σ = 1.2
f(u, p, t) = μ * u
g(u, p, t) = σ * u
prob = SDEProblem(f, g, 1.0, (0.0, Inf))
```

Now we build the noise process by building the integrator and sending that
integrator to the `NoiseApproximation` constructor:

```julia
integrator = init(prob, SRIW1())
W = NoiseApproximation(integrator)
```

We can use this noise process like any other noise process. For example, we
can now build a geometric Brownian motion whose noise process is colored noise
that itself is a geometric Brownian motion:

```julia
prob = SDEProblem(f, g, 1.0, (0.0, Inf), noise = W)
```

The possibilities are endless.
"""
mutable struct NoiseApproximation{T, N, Tt, T2, T3, S1, S2, ZType, inplace} <:
               AbstractNoiseProcess{T, N, Vector{T2}, inplace}
    t::Vector{Tt}
    u::Vector{T2}
    W::Vector{T2}
    Z::ZType
    curt::Tt
    curW::T2
    curZ::T3
    dt::Tt
    dW::T2
    dZ::T3
    source1::S1
    source2::S2
    reset::Bool
end

function NoiseApproximation(source1::DEIntegrator,
        source2::Union{DEIntegrator, Nothing} = nothing;
        reset = true)
    _source1 = deepcopy(source1)
    _source2 = deepcopy(source2)
    if _source2 == nothing
        Z = nothing
        curZ = nothing
        dZ = nothing
    else
        Z = _source2.sol.u
        curZ = copy(_source2.u)
        dZ = copy(_source2.u)
        _source2.opts.advance_to_tstop = true
    end
    val = copy(_source1.u)
    t = _source1.sol.t
    W = _source1.sol.u
    curW = copy(_source1.u)
    dW = copy(_source1.u)
    dt = _source1.dt
    curt = _source1.t
    _source1.opts.advance_to_tstop = true
    NoiseApproximation{typeof(val), ndims(val), typeof(curt), typeof(curW), typeof(curZ),
        typeof(_source1), typeof(_source2), typeof(Z),
        isinplace(_source1.sol.prob)}(t, W, W, Z, curt, curW, curZ, dt, dW,
        dZ, _source1, _source2, reset)
end

(W::NoiseApproximation)(t) = interpolate!(W, t)
(W::NoiseApproximation)(u, p, t) = interpolate!(W, t)
(W::NoiseApproximation)(out1, out2, u, p, t) = interpolate!(out1, out2, W, t)

"""
A `VirtualBrownianTree` builds the noise process starting from an initial
time `t0`, the first value of the process `W0`, and (optionally) the first
value `Z0` for an auxiliary pseudo-process. The constructor is given as

```julia
VirtualBrownianTree(t0,
    W0,
    Z0 = nothing,
    dist = WHITE_NOISE_DIST,
    bridge = VBT_BRIDGE;
    kwargs...)
```

where `dist` specifies the distribution that is used to generate the end
point(s) `Wend` (`Zend`) of the noise process for the final time `tend`.
`bridge` denotes the distribution of the employed Brownian bridge.  Per
default `tend` is fixed to `t0+1` but can be changed by passing a custom `tend`
as a keyword argument. The following keyword arguments are available:

  - `tend` is the end time of the noise process.
  - `Wend` is the end value of the noise process.
  - `Zend` is the end value of the pseudo-noise process.
  - `atol` represents the absolute tolerance determining when the recursion is
    terminated.
  - `tree_depth` allows one to store a cache of seeds, noise values, and times
    to speed up the simulation by reducing the recursion steps.
  - `search_depth` maximal search depth for the tree if `atol` is not reached.
  - `rng` the splittable PRNG used for generating the random numbers.
    Default: `Threefry4x()` from the Random123 package.

## VirtualBrownianTree Example

In this example, we define a mulitidimensional Brownian process based on a
`VirtualBrownianTree` with a minimal `tree_depth=0` such that memory consumption
is minimized.

```julia
W0 = zeros(10)
W = VirtualBrownianTree(0.0, W0; tree_depth = 0)

prob = NoiseProblem(W, (0.0, 1.0))
sol = solve(prob; dt = 1 / 10)
```

Using a look-up cache by increasing `tree_depth` can significantly reduce the
runtime. Thus, the `VirtualBrownianTree` allows for trading off speed for memory
in a simple manner.
"""
mutable struct VirtualBrownianTree{
    T, N, F, F2, Tt, T2, T3, T2tmp, T3tmp, seedType, tolType,
    RNGType, inplace} <:
               AbstractNoiseProcess{T, N, Vector{T2}, inplace}
    dist::F
    bridge::F2
    t::Vector{Tt} # tstart::Tt to tend::Tt
    u::Vector{T2}
    W::Vector{T2}
    Z::Vector{T3}
    curt::Tt
    curW::T2
    curZ::T3
    dt::Tt
    dW::T2
    dZ::T3
    W0tmp::T2tmp
    W1tmp::T2tmp
    Z0tmp::T3tmp
    Z1tmp::T3tmp
    seeds::Vector{seedType}
    atol::tolType
    rng::RNGType
    tree_depth::Int
    search_depth::Int
    step_setup::Bool
    reset::Bool
end

function VirtualBrownianTree{iip}(t0, W0, Z0, dist, bridge;
        tend = nothing, Wend = nothing, Zend = nothing,
        atol = 1e-10, tree_depth::Int = 4,
        search_depth = nothing,
        rng = RandomNumbers.Random123.Threefry4x(),
        reset = true) where {iip}
    if search_depth == nothing
        if atol < 1e-10
            search_depth = 50 # maximum search depth
        else
            search_depth = floor(Int, (-log2(atol) + 2))
        end
    end

    (tree_depth >= search_depth) &&
        error("search depth of the VBT must be greater than the depth of the cached tree")

    curt = t0
    dt = t0
    curW = copy(W0)
    dW = copy(W0)
    if Z0 == nothing
        curZ = nothing
        dZ = nothing
    else
        curZ = copy(Z0)
        dZ = copy(Z0)
    end

    # assign final time and state
    if tend == nothing
        tend = t0 + one(t0)
    end

    if Wend == nothing
        Wend = curW + dist(dW, nothing, tend - t0, nothing, nothing, nothing, rng)
    end

    if Zend == nothing && Z0 !== nothing
        Zend = curW + bridge(dZ, nothing, tend - t0, nothing, nothing, nothing, rng)
    end

    t, W, Z, seeds = create_VBT_cache(
        bridge, t0, W0, Z0, tend, Wend, Zend, rng, tree_depth,
        search_depth)

    if iip
        W0tmp, W1tmp = copy(W0), copy(Wend)
        if Z0 !== nothing
            Z0tmp, Z1tmp = copy(Z0), copy(Zend)
        else
            Z0tmp, Z1tmp = nothing, nothing
        end
    else
        W0tmp, W1tmp, Z0tmp, Z1tmp = nothing, nothing, nothing, nothing
    end

    VirtualBrownianTree{
        typeof(W0), ndims(W0),
        typeof(dist), typeof(bridge), typeof(curt), typeof(curW),
        typeof(curZ),
        typeof(W0tmp), typeof(Z0tmp), typeof(seeds[1]), typeof(atol),
        typeof(rng),
        iip}(dist, bridge, t, W, W, Z, curt, curW, curZ, dt, dW, dZ, W0tmp,
        W1tmp, Z0tmp, Z1tmp, seeds, atol, rng, tree_depth,
        search_depth, true, reset)
end

(W::VirtualBrownianTree)(t) = interpolate!(W, nothing, nothing, t)
(W::VirtualBrownianTree)(u, p, t) = interpolate!(W, u, p, t)
(W::VirtualBrownianTree)(out1, out2, u, p, t) = interpolate!(out1, out2, W, u, p, t)

function VirtualBrownianTree(t0, W0, Z0 = nothing, dist = WHITE_NOISE_DIST,
        bridge = VBT_BRIDGE; kwargs...)
    VirtualBrownianTree{false}(t0, W0, Z0, dist, bridge; kwargs...)
end

function VirtualBrownianTree!(t0, W0, Z0 = nothing, dist = INPLACE_WHITE_NOISE_DIST,
        bridge = INPLACE_VBT_BRIDGE; kwargs...)
    VirtualBrownianTree{true}(t0, W0, Z0, dist, bridge; kwargs...)
end

"""
```julia
BoxWedgeTail{T, N, Tt, TA, T2, T3, ZType, F, F2, inplace, RNGType, tolType, spacingType,
    jpdfType, boxType, wedgeType, tailType, distBWTType, distΠType} <:
AbstractNoiseProcess{T, N, Vector{T2}, inplace}
```

The method for random generation of stochastic area integrals due to Gaines and Lyons. The method is
based on Marsaglia's "rectangle-wedge-tail" approach for two dimensions.

3 different groupings for the boxes are implemented.

  - box_grouping = :Columns (full, i.e., as large as possible, columns on a square spanned by dr and da)
  - box_grouping = :none (no grouping)
  - box_grouping = :MinEntropy (default, grouping that achieves a smaller entropy than the column wise grouping and thus allows for slightly faster sampling -- but has a slightly larger number of groups)

The sampling is based on the Distributions.jl package, i.e., to sample from one of the many distributions,
a uni-/bi-variate distribution from Distributions.jl is constructed, and then rand(..) is used.

## Constructor

```julia
BoxWedgeTail{iip}(t0, W0, Z0, dist, bridge;
    rtol = 1e-8, nr = 4, na = 4, nz = 10,
    box_grouping = :MinEntropy,
    sqeezing = true,
    save_everystep = true,
    rng = Xorshifts.Xoroshiro128Plus(rand(UInt64)),
    reset = true, reseed = true) where {iip}
```
"""
mutable struct BoxWedgeTail{T, N, Tt, TA, T2, T3, ZType, F, F2, inplace, RNGType, tolType,
    spacingType, jpdfType, boxType, wedgeType, tailType,
    distBWTType, distΠType} <:
               AbstractNoiseProcess{T, N, Vector{T2}, inplace}
    dist::F
    bridge::F2
    t::Vector{Tt}
    A::Vector{TA}
    u::Vector{T2} # Aliased pointer to W for the AbstractVectorOfArray interface
    W::Vector{T2}
    Z::ZType
    curt::Tt
    curA::TA
    curW::T2
    curZ::T3
    dt::Tt
    dW::T2
    dZ::T3
    dWtilde::T2
    dZtilde::T3
    dWtmp::T2
    dZtmp::T3
    save_everystep::Bool
    iter::Int
    rng::RNGType
    reset::Bool
    reseed::Bool
    rtol::tolType
    Δr::spacingType
    Δa::spacingType
    Δz::spacingType
    rM::spacingType
    aM::spacingType
    jpdf::jpdfType
    boxes::boxType
    wedges::wedgeType
    sqeezing::Bool
    tails::tailType
    distBWT::distBWTType
    distΠ::distΠType
end

function BoxWedgeTail{iip}(t0, W0, Z0, dist, bridge;
        rtol = 1e-8, nr = 4, na = 4, nz = 10,
        box_grouping = :MinEntropy,
        sqeezing = true,
        save_everystep = true,
        rng = Xorshifts.Xoroshiro128Plus(rand(UInt64)),
        reset = true, reseed = true) where {iip}
    if Z0 === nothing
        Z = nothing
        curZ = nothing
        dZ = nothing
        dZtilde = nothing
        dZtmp = nothing
    else
        Z = [copy(Z0)]
        curZ = copy(Z0)
        dZ = copy(Z0)
        dZtilde = copy(Z0)
        dZtmp = copy(Z0)
    end
    W = [copy(W0)]
    A0 = zero(eltype(W0))
    N = length((size(W0)..., length(W)))

    N != 2 &&
        error("The BoxWedgeTail algorithms can only be used for 2d Brownian processes.")

    jpdf = (r, a) -> joint_density_function(r, a, rtol)

    # grid in a
    aM = 4 * one(eltype(W0))
    Δa = convert(typeof(aM), 2)^(-na)

    # grid in r
    rM = 4 * one(eltype(W0))
    Δr = convert(typeof(rM), 2)^(-nr)

    # smallest z value
    Δz = convert(typeof(rM), 2)^(-nz)

    # generate boxes
    if box_grouping == :MinEntropy
        box, probability, offset = generate_boxes2(jpdf, Δr, Δa, Δz, one(Δr), one(Δa),
            one(Δz) / 64, rM, aM)
        dist_box = Distributions.Categorical(probability / sum(probability))
        boxes = BoxGeneration2{typeof(box), typeof(probability), typeof(offset),
            typeof(dist_box)}(box,
            probability,
            offset,
            dist_box)
    elseif box_grouping == :Columns
        box, probability, offset = generate_boxes1(jpdf, Δr, Δa, Δz, rM, aM)
        val1 = sum(probability)
        dist_box = Distributions.Categorical(probability / val1)
        boxes = BoxGeneration1{typeof(box), typeof(probability), typeof(offset),
            typeof(dist_box)}(box,
            probability,
            offset,
            dist_box)
    elseif box_grouping == :none
        box, probability, offset = generate_boxes3(jpdf, Δr, Δa, Δz, rM, aM)
        val1 = sum(probability)
        dist_box = Distributions.Categorical(probability / val1)
        boxes = BoxGeneration3{typeof(box), typeof(probability), typeof(offset),
            typeof(dist_box)}(box,
            probability,
            offset,
            dist_box)
    else
        error("Available options for box grouping are :MinEntropy, :Columns, and :none.")
    end

    # generate wedges
    box, probability, offset = generate_wedges(jpdf, Δr, Δa, Δz, rM, aM, offset, sqeezing)
    dist_box = Distributions.Categorical(probability / sum(probability))
    wedges = Wedges{typeof(box), typeof(probability), typeof(offset),
        typeof(dist_box)}(box,
        probability,
        offset,
        dist_box)

    # set up tail approximation
    tails = TailApproxs(rM, aM)

    # distribution to decide if sample should be drawn from boxes, wedges or tail
    if nr == 4 && na == 4 && nz == 10
        distBWT = Distributions.Categorical([
            0.9118576049804688,
            0.08546645498798722,
            0.002675940031544033
        ])
    else
        error("Cubature not implemented but needed for a different discretization. Please report this error.")
    end

    # distribution between 0 and 2pi to convert r to dWs
    distΠ = Distributions.Uniform(zero(rM), convert(typeof(rM), 2 * pi))

    BoxWedgeTail{eltype(eltype(W0)), N, typeof(t0), typeof(A0), typeof(W0), typeof(dZ),
        typeof(Z),
        typeof(dist), typeof(bridge), iip, typeof(rng), typeof(Δr), typeof(rtol),
        typeof(jpdf), typeof(boxes), typeof(wedges), typeof(tails),
        typeof(distBWT), typeof(distΠ)}(dist, bridge, [t0], [A0], W, W, Z, t0, A0,
        copy(W0), curZ, t0, copy(W0), dZ, copy(W0),
        dZtilde, copy(W0), dZtmp,
        save_everystep, 0, rng, reset, reseed,
        rtol, Δr, Δa, Δz, rM, aM, jpdf, boxes,
        wedges,
        sqeezing, tails, distBWT, distΠ)
end

(W::BoxWedgeTail)(t) = interpolate!(W, nothing, nothing, t)
(W::BoxWedgeTail)(u, p, t) = interpolate!(W, u, p, t)
(W::BoxWedgeTail)(out1, out2, u, p, t) = interpolate!(out1, out2, W, u, p, t)

function BoxWedgeTail(t0, W0, Z0 = nothing, dist = WHITE_NOISE_DIST,
        bridge = WHITE_NOISE_BRIDGE; kwargs...)
    BoxWedgeTail{false}(t0, W0, Z0, dist, bridge; kwargs...)
end

function BoxWedgeTail!(t0, W0, Z0 = nothing, dist = INPLACE_WHITE_NOISE_DIST,
        bridge = INPLACE_WHITE_NOISE_BRIDGE; kwargs...)
    BoxWedgeTail{true}(t0, W0, Z0, dist, bridge; kwargs...)
end
