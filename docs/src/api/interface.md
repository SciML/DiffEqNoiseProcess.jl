# Noise Process Interface

This page documents the interface functions for working with noise processes.

## Process Contract

All noise objects in this package are subtypes of
[`SciMLBase.AbstractNoiseProcess`](https://docs.sciml.ai/SciMLBase/stable/interfaces/Solutions/#SciMLBase.AbstractNoiseProcess).
The four type parameters are the scalar element type, array rank, saved-value
array type, and the `isinplace` mutation convention. A concrete process must
also behave like an `AbstractDiffEqArray`: saved values are indexed through the
usual array interface, and callable processes evaluate `(W, Z)` at a time.

The callable forms are:

```julia
W(t)                    # returns (primary_value, auxiliary_value)
W(u, p, t)              # state- and parameter-aware form
W(out1, out2, u, p, t)  # in-place form when isinplace(W) is true
```

When no auxiliary process exists, the second returned value is `nothing`. The
in-place form writes into `out1` and, when applicable, `out2`; it must not
replace those caller-owned buffers.

## Solver Lifecycle

SDE and noise solvers use the exported lifecycle functions in this order:

1. `setup_next_step!(W, u, p)` prepares the pending increment for `W.dt`.
2. `calculate_step!(W, dt, u, p)` recomputes a pending increment without
   committing it.
3. `accept_step!(W, dt, u, p, setup_next)` commits the pending increment once
   and optionally prepares the next one.
4. `reject_step!(W, dtnew, u, p)` replaces a rejected increment while preserving
   its conditional distribution. A process that cannot support rejection must
   raise an explicit error; it must not silently reuse the old increment.
5. `save_noise!(W)` records the current point for history-backed processes.

The `u` and `p` arguments are `nothing` for state-independent noise. A process
with state-dependent callbacks must forward them unchanged. Implementations
must keep primary and auxiliary state synchronized. `accept_step!`,
`reject_step!`, `setup_next_step!`, and `save_noise!` return `nothing`;
`calculate_step!` may additionally return the selected step width for
grid-backed processes, but callers should use the process state.

The step functions are developer-facing extension points for solver packages;
the process's stack, cache, and scratch-buffer fields are not part of this
contract. Use the lifecycle functions rather than mutating those fields.

## Step Management

```@docs
accept_step!
reject_step!
calculate_step!
setup_next_step!
save_noise!
```

## Configuration

### Rejection Sampling with Memory (RSWM)

```@docs
RSWM
adaptive_alg
```

## Developer Callback Helpers

The following callback helpers are documented for package developers who need to
understand or compose the built-in process implementations. They are not the
stable user-facing process API and are not exported. New solver code should use
the lifecycle functions above.

### Distribution Functions

```@docs
DiffEqNoiseProcess.WHITE_NOISE_DIST
DiffEqNoiseProcess.WHITE_NOISE_BRIDGE
DiffEqNoiseProcess.VBT_BRIDGE
DiffEqNoiseProcess.INPLACE_WHITE_NOISE_DIST
DiffEqNoiseProcess.INPLACE_WHITE_NOISE_BRIDGE
DiffEqNoiseProcess.INPLACE_VBT_BRIDGE
DiffEqNoiseProcess.REAL_WHITE_NOISE_DIST
DiffEqNoiseProcess.REAL_WHITE_NOISE_BRIDGE
DiffEqNoiseProcess.REAL_INPLACE_WHITE_NOISE_DIST
DiffEqNoiseProcess.REAL_INPLACE_WHITE_NOISE_BRIDGE
```

### Random Number Generation

```@docs
DiffEqNoiseProcess.wiener_randn
DiffEqNoiseProcess.wiener_randn!
```

### Ornstein-Uhlenbeck Specific

```@docs
DiffEqNoiseProcess.OrnsteinUhlenbeck
DiffEqNoiseProcess.OrnsteinUhlenbeck!
DiffEqNoiseProcess.ou_bridge
DiffEqNoiseProcess.ou_bridge!
```

### Geometric Brownian Motion Specific

```@docs
DiffEqNoiseProcess.GeometricBrownianMotion
DiffEqNoiseProcess.GeometricBrownianMotion!
DiffEqNoiseProcess.gbm_bridge
DiffEqNoiseProcess.gbm_bridge!
```

### Compound Poisson Specific

```@docs
DiffEqNoiseProcess.cpp_bridge
DiffEqNoiseProcess.cpp_bridge!
```
