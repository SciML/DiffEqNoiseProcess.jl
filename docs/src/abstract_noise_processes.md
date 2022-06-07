# Abstract Noise Processes

In addition to the `NoiseProcess` type, more general `AbstractNoiseProcess`es
are defined. The `NoiseGrid` allows you to define a noise process from a set
of pre-calculated points (the "normal" way). The `NoiseApproximation` allows
you to define a new noise process as the solution to some stochastic differential
equation. While these methods are only approximate, they are more general and
allow the user to easily define their own colored noise to use in simulations.

The `NoiseWrapper` allows one to wrap a `NoiseProcess` from a previous simulation
to re-use it in a new simulation in a way that follows the same stochastic
trajectory (even if different points are hit, for example solving with a
smaller `dt`) in a distributionally-exact manner. It is demonstrated how the
`NoiseWrapper` can be used to wrap the `NoiseProcess` of one SDE/RODE solution
in order to re-use the same noise process in another simulation.

The `VirtualBrownianTree` allows one to trade speed for O(1) memory usage.
Instead of storing Brownian motion increments, the `VirtualBrownianTree` samples
recursively from the midpoint `tmid` of Brownian bridges, using a splittable PRNG.
The recursion terminates when the query time agrees within some tolerance
with `tmid` or when the maximum depth of the tree is reached.

Lastly, the `NoiseFunction` allows you to use any function of time as the
noise process. Together, this functionality allows you to define any colored
noise process and use this efficiently and accurately in your simulations.

## The Standard `AbstractNoiseProcess`

```@docs
NoiseProcess
```

## Alternative `AbstractNoiseProcess` Types

In addition to the mathematically-defined noise processes above, there exists
more generic functionality for building noise processes from other noise processes,
from arbitrary functions, from arrays, and from approximations of stochastic
differential equations.

```@docs
NoiseWrapper
NoiseFunction
NoiseGrid
NoiseApproximation
VirtualBrownianTree
SimpleNoiseProcess
BoxWedgeTail
```