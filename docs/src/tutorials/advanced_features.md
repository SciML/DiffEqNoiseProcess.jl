# Advanced Features

This tutorial covers the advanced functionality in DiffEqNoiseProcess.jl, including custom noise processes, wrappers, and memory-efficient alternatives.

## Custom Noise from Functions

### NoiseFunction

Create noise processes from arbitrary functions of time:

```@example advanced
using DiffEqNoiseProcess, SciMLBase
using Random

# Define a deterministic "noise" function
# NoiseFunction expects signature f(u, p, t)
function my_noise_func(u, p, t)
    return sin(2π * t) + 0.1 * cos(10π * t)
end

# Create noise process from function (t0, function)
noise_func = NoiseFunction(0.0, my_noise_func)

# Test the function
prob = NoiseProblem(noise_func, (0.0, 2.0))
sol = solve(prob; dt = 0.01)

println("NoiseFunction example completed successfully!")
# Access values via callable interface: sol(t) returns (W, Z) tuple
println("Value at t=0.5: $(sol(0.5)[1])")
println("Value at t=1.0: $(sol(1.0)[1])")
```

### NoiseTransport

Transport a random variable through a time-dependent function:

```@example advanced
# Transport function: signature is f(u, p, t, rv) where rv is the random variable
# This creates noise W(t) = ξ * exp(-0.5t) * cos(t) where ξ ~ N(0,1)
function transport_func(u, p, t, ξ)
    return ξ * exp(-0.5 * t) * cos(t)
end

# Create transported noise: NoiseTransport(t0, f, RV)
# RV is a function that generates random values (like randn)
noise_transport = NoiseTransport(0.0, transport_func, randn)

prob_transport = NoiseProblem(noise_transport, (0.0, 3π))
sol_transport = solve(prob_transport; dt = 0.01)

# NoiseTransport computes values on-the-fly; access via callable interface
# Returns (W_value, Z_value) tuple, use [1] to get W
println("Transported at t=0: $(sol_transport(0.0)[1])")
println("Transported at t=π: $(sol_transport(π)[1])")
```

## Noise from Data

### NoiseGrid

Create noise processes from pre-computed data points:

```@example advanced
# Generate some data (e.g., from experimental measurements)
t_data = collect(0.0:0.1:1.0)
noise_data = cumsum(randn(length(t_data)) * sqrt(0.1))  # Approximate Brownian motion

# Create noise process from grid
noise_grid = NoiseGrid(t_data, noise_data)

# Use it in a problem
prob_grid = NoiseProblem(noise_grid, (0.0, 1.0))
sol_grid = solve(prob_grid; dt = 0.05)

println("Grid noise interpolated values:")
for t in [0.0, 0.25, 0.5, 0.75, 1.0]
    println("  t=$t: $(sol_grid(t))")
end
```

## Memory-Efficient Alternatives

### VirtualBrownianTree

For memory-constrained applications, use a virtual Brownian tree that generates values on-demand:

```@example advanced
# Create a virtual Brownian tree
# VirtualBrownianTree(t0, W0; tree_depth, atol, ...)
# tree_depth controls the cache size for speed/memory tradeoff
vbt = VirtualBrownianTree(0.0, 0.0; tree_depth = 5, atol = 1e-6)

# The tree generates Brownian motion values on demand without storing them
prob_vbt = NoiseProblem(vbt, (0.0, 1.0))
sol_vbt = solve(prob_vbt; dt = 0.01)

println("VBT memory usage is O(tree_depth) regardless of path length")
println("Final VBT value: $(sol_vbt.u[end])")
```

## Noise Process Wrappers

### NoiseWrapper

Reuse noise from previous simulations in a distributionally-exact manner:

```@example advanced
# First, generate a reference noise process
W_ref = WienerProcess(0.0, 0.0, 1.0)
prob_ref = NoiseProblem(W_ref, (0.0, 1.0))
sol_ref = solve(prob_ref; dt = 0.1)  # Coarse timestep

println("Reference solution has $(length(sol_ref.t)) points")

# Now wrap it to use at different timesteps
# Note: solve(NoiseProblem(...)) returns the noise process directly
wrapped_noise = NoiseWrapper(sol_ref)

# Use the wrapped noise with finer timestep
prob_wrapped = NoiseProblem(wrapped_noise, (0.0, 1.0))
sol_wrapped = solve(prob_wrapped; dt = 0.01)  # Fine timestep

println("Wrapped solution has $(length(sol_wrapped.t)) points")
println("Both solutions follow the same stochastic trajectory")

# Verify they agree at common timepoints
# Note: WienerProcess(t0, W0, Z0) returns tuples (W, Z), so we access [1] for W
for t in [0.0, 0.5, 1.0]
    ref_val = sol_ref(t)[1]
    wrapped_val = sol_wrapped(t)[1]
    println("At t=$t: ref=$(round(ref_val, digits=6)), wrapped=$(round(wrapped_val, digits=6))")
end
```

## Approximation Methods

### NoiseApproximation

Approximate colored noise as solutions to SDEs. `NoiseApproximation` takes a
`DEIntegrator` from an SDE solve, allowing you to use SDE solutions as noise
processes. This is useful for creating correlated or colored noise.

```@example advanced
# NoiseApproximation requires StochasticDiffEq for the SDE integrator
# Here we show a conceptual example - in practice you would do:
#
# using StochasticDiffEq
# f(u, p, t) = -u          # Mean reversion (drift)
# g(u, p, t) = 0.5         # Constant diffusion
# sde_prob = SDEProblem(f, g, 1.0, (0.0, Inf))
# integrator = init(sde_prob, SRIW1())
# noise_approx = NoiseApproximation(integrator)
#
# Then use noise_approx as a noise process in another problem

println("NoiseApproximation requires a DEIntegrator from StochasticDiffEq")
println("See the API documentation for full usage examples")
```

## Process Configuration

### Process Reset and Reuse

```@example advanced
# Create a process and use it
W_reusable = WienerProcess(0.0, 0.0, 1.0)
prob1 = NoiseProblem(W_reusable, (0.0, 1.0))
sol1 = solve(prob1; dt = 0.1)

println("First use - final value: $(sol1.u[end])")

# Reinitialize and reuse the same process
reinit!(W_reusable, 0.1)  # Reinitialize with dt=0.1
prob2 = NoiseProblem(W_reusable, (0.0, 1.0))
sol2 = solve(prob2; dt = 0.1)

println("After reinit - final value: $(sol2.u[end])")
println("Values are different due to reinitialization")
```

## Performance Considerations

Choose the appropriate noise process type based on your application needs:

1. **In-place versions**: Use in-place versions (Process!) for large systems to reduce allocations
2. **VirtualBrownianTree**: Use for memory-constrained applications where O(1) memory is critical
3. **NoiseGrid**: Use for pre-computed data or when you have experimental measurements
4. **NoiseWrapper**: Use to reuse expensive computations across multiple simulations with different timesteps