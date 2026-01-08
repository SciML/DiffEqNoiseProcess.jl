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
println("Solution contains $(length(sol.t)) time points")
```

### NoiseTransport

Transport a random variable through a time-dependent function:

```@example advanced
using Distributions

# Random initial value
Random.seed!(456)
ξ = rand(Normal(0, 1))  # Standard normal random variable

# Transport function: moves and scales the random variable over time
function transport_func(ξ_val, t)
    return ξ_val * exp(-0.5 * t) * cos(t)
end

# Create transported noise
noise_transport = NoiseTransport(ξ, transport_func, 0.0)

prob_transport = NoiseProblem(noise_transport, (0.0, 3π))
sol_transport = solve(prob_transport; dt = 0.01)

println("Initial random value: $ξ")
println("Transported at t=0: $(sol_transport.u[1])")
println("Transported at t=π: $(sol_transport(π))")
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
using Random

# Create a virtual Brownian tree
# Parameters: t0, T (end time), tree_depth, tolerance, rng
rng = MersenneTwister(789)
vbt = VirtualBrownianTree(0.0, 1.0, 10, 1e-6, rng)

# The tree generates Brownian motion values on demand without storing them
prob_vbt = NoiseProblem(vbt, (0.0, 1.0))
sol_vbt = solve(prob_vbt; dt = 0.01)

println("VBT memory usage is O(1) regardless of path length")
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
wrapped_noise = NoiseWrapper(sol_ref.W)

# Use the wrapped noise with finer timestep
prob_wrapped = NoiseProblem(wrapped_noise, (0.0, 1.0))
sol_wrapped = solve(prob_wrapped; dt = 0.01)  # Fine timestep

println("Wrapped solution has $(length(sol_wrapped.t)) points")
println("Both solutions follow the same stochastic trajectory")

# Verify they agree at common timepoints
for t in [0.0, 0.5, 1.0]
    ref_val = sol_ref(t)
    wrapped_val = sol_wrapped(t)
    println("At t=$t: ref=$(round(ref_val, digits=6)), wrapped=$(round(wrapped_val, digits=6))")
end
```

## Approximation Methods

### NoiseApproximation

Approximate colored noise as solutions to SDEs:

```@example advanced
# Define an SDE that generates colored noise
# For example, dX = -X*dt + σ*dW (Ornstein-Uhlenbeck-like)
function colored_drift(u, p, t)
    return -u  # Mean reversion
end

function colored_diffusion(u, p, t)
    return 0.5  # Constant diffusion
end

# Create approximation
noise_approx = NoiseApproximation(colored_drift, colored_diffusion, 0.0, 0.0, 1.0)

prob_approx = NoiseProblem(noise_approx, (0.0, 2.0))
sol_approx = solve(prob_approx; dt = 0.01)

println("Colored noise approximation:")
println("  Initial value: $(sol_approx.u[1])")
println("  Final value: $(sol_approx.u[end])")
```

## Process Configuration

### Process Reset and Reuse

```@example advanced
# Create a process and use it
W_reusable = WienerProcess(0.0, 0.0, 1.0)
prob1 = NoiseProblem(W_reusable, (0.0, 1.0))
sol1 = solve(prob1; dt = 0.1)

println("First use - final value: $(sol1.u[end])")

# Reset and reuse the same process type
reset!(W_reusable)  # Reset internal state
prob2 = NoiseProblem(W_reusable, (0.0, 1.0))
sol2 = solve(prob2; dt = 0.1)

println("After reset - final value: $(sol2.u[end])")
println("Values are different due to reset")
```

## Performance Considerations

Choose the appropriate noise process type based on your application needs:

1. **In-place versions**: Use in-place versions (Process!) for large systems to reduce allocations
2. **VirtualBrownianTree**: Use for memory-constrained applications where O(1) memory is critical
3. **NoiseGrid**: Use for pre-computed data or when you have experimental measurements
4. **NoiseWrapper**: Use to reuse expensive computations across multiple simulations with different timesteps