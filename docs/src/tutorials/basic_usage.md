# Basic Usage Tutorial

This tutorial covers the fundamental usage patterns of DiffEqNoiseProcess.jl, including creating noise processes, using them in problems, and direct simulation.

## Creating and Using Noise Processes

### Basic Wiener Process

A Wiener process (Brownian motion) is the most fundamental noise process:

```@example basic
using DiffEqNoiseProcess, SciMLBase

# Create a standard Wiener process
# Parameters: t0 (start time), W0 (initial value), Z0 (end value)
W = WienerProcess(0.0, 0.0, 1.0)

# You can specify additional options
W_custom = WienerProcess(0.0, 0.0, 1.0; 
                        save_everystep = true,
                        rswm = RSWM())
```

### Real-Valued Wiener Process

For scalar problems, use a real-valued Wiener process:

```@example basic
# Real-valued Wiener process  
W_real = RealWienerProcess(0.0, 0.0, 1.0)
```

### Geometric Brownian Motion

A geometric Brownian motion process is useful for financial modeling:

```@example basic
# Parameters: μ (drift), σ (volatility), t0, W0, Z0
μ = 0.05  # 5% drift
σ = 0.2   # 20% volatility
GBM = GeometricBrownianMotionProcess(μ, σ, 0.0, 1.0, 1.0)
```

## Using Noise in SDE Problems

Noise processes can be passed directly to SDE problems:

```@example basic
# Define a simple SDE: du = u*dt + 0.1*u*dW
function f(u, p, t)
    return u
end

function g(u, p, t)
    return 0.1 * u
end

u0 = 1.0
tspan = (0.0, 1.0)

# Create the problem with custom noise
prob = SDEProblem(f, g, u0, tspan, noise = W)
```

## Direct Simulation

You can simulate noise processes directly without solving an SDE:

```@example basic
# Create a noise problem
prob = NoiseProblem(W, (0.0, 1.0))

# Solve it with a specified timestep
sol = solve(prob; dt = 0.01)

# Access the noise values
println("Final noise value: ", sol.u[end])
println("Number of steps: ", length(sol.t))
```

### Ensemble Simulations

Generate multiple realizations of the same noise process:

```@example basic
# Create ensemble problem
enprob = EnsembleProblem(prob)

# Generate 100 trajectories
ensol = solve(enprob; dt = 0.01, trajectories = 100)

println("Generated $(length(ensol.u)) noise trajectories")
```

## Direct Interface for Advanced Usage

For advanced users, noise processes support a direct stepping interface:

```@example basic
# Create a fresh Wiener process
W_direct = WienerProcess(0.0, 0.0, 1.0)

# Set the timestep
dt = 0.1
W_direct.dt = dt

# For state-dependent noise (not applicable here, but required by interface)
u = nothing
p = nothing

# Calculate and accept steps
calculate_step!(W_direct, dt, u, p)
for i in 1:5
    accept_step!(W_direct, dt, u, p)
    println("Step $i: W = $(W_direct.curW)")
end
```

This direct interface is primarily used internally by the SDE solvers but can be useful for custom implementations.