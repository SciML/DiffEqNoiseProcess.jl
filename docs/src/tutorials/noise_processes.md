# Noise Process Types

This tutorial demonstrates the various types of noise processes available in DiffEqNoiseProcess.jl and their specific use cases.

## Standard Brownian Motion

### Wiener Process

The standard Wiener process is the foundation of stochastic calculus:

```@example noise_types
using DiffEqNoiseProcess, SciMLBase, Random, Distributions
using Statistics

# Set seed for reproducibility
Random.seed!(123)

# Standard Wiener process
W = WienerProcess(0.0, 0.0, 1.0)

# Simulate it
prob = NoiseProblem(W, (0.0, 2.0))
sol = solve(prob; dt = 0.01)

println("Final value: $(sol.u[end])")
println("Path length: $(length(sol.t))")
```

### Correlated Wiener Process

For multi-dimensional systems with correlated noise:

```@example noise_types
# 2D correlated Wiener process with correlation matrix
Γ = [1.0 0.5; 0.5 1.0]  # Correlation matrix
W0 = [0.0, 0.0]         # Initial values
Z0 = [1.0, 1.0]         # Final values

corr_W = CorrelatedWienerProcess(Γ, 0.0, W0, Z0)

# Simulate
prob_corr = NoiseProblem(corr_W, (0.0, 1.0))
sol_corr = solve(prob_corr; dt = 0.01)

println("Final 2D noise: $(sol_corr.u[end])")
```

## Mean-Reverting Processes

### Ornstein-Uhlenbeck Process

Models mean-reverting behavior, which is commonly used in interest rate models:

```@example noise_types
# Parameters: Θ (mean reversion rate), μ (long-term mean), σ (volatility)
Θ = 2.0   # Fast mean reversion
μ = 0.5   # Mean level
σ = 0.3   # Volatility

OU = OrnsteinUhlenbeckProcess(Θ, μ, σ, 0.0, 0.0, 1.0)

# Simulate
prob_ou = NoiseProblem(OU, (0.0, 5.0))
sol_ou = solve(prob_ou; dt = 0.01)

println("OU process starts at: $(sol_ou.u[1])")
println("OU process ends at: $(sol_ou.u[end])")
println("Long-term mean μ: $μ")
```

## Jump Processes

### Compound Poisson Process

Models discontinuous jumps at random times:

```@example noise_types
# Parameters: λ (jump rate), jump distribution
λ = 2.0  # 2 jumps per unit time on average

# Jump sizes are normally distributed N(0, 0.1²)
jump_dist = Normal(0.0, 0.1)

poisson_proc = CompoundPoissonProcess(λ, 0.0, 0.0)

# Simulate
prob_poisson = NoiseProblem(poisson_proc, (0.0, 2.0))
sol_poisson = solve(prob_poisson; dt = 0.01)

println("Final Poisson process value: $(sol_poisson.u[end])")

# Count number of jumps (non-zero increments)
jumps = sum(abs.(diff(sol_poisson.u)) .> 1e-10)
println("Number of jumps: $jumps")
```

## Financial Models

### Geometric Brownian Motion

The classic Black-Scholes model for stock prices:

```@example noise_types
# Stock parameters
μ_stock = 0.08   # 8% expected return
σ_stock = 0.25   # 25% volatility
S0 = 100.0       # Initial stock price

GBM = GeometricBrownianMotionProcess(μ_stock, σ_stock, 0.0, S0, S0)

# Simulate stock price over 1 year
prob_gbm = NoiseProblem(GBM, (0.0, 1.0))
sol_gbm = solve(prob_gbm; dt = 1/252)  # Daily steps

println("Initial stock price: $(sol_gbm.u[1])")
println("Final stock price: $(sol_gbm.u[end])")
println("Return: $((sol_gbm.u[end] - sol_gbm.u[1])/sol_gbm.u[1] * 100)%")
```

## Bridge Processes

Bridge processes are conditioned to pass through specific values:

### Brownian Bridge

A Wiener process conditioned to start and end at specific values:

```@example noise_types
# Bridge from 0 at t=0 to 1 at t=1
bridge = BrownianBridge(0.0, 1.0, 0.0, 1.0)

prob_bridge = NoiseProblem(bridge, (0.0, 1.0))
sol_bridge = solve(prob_bridge; dt = 0.01)

println("Bridge starts at: $(sol_bridge.u[1])")
println("Bridge ends at: $(sol_bridge.u[end])")
println("Should be exactly 0 and 1 respectively")
```

### Geometric Brownian Bridge

A geometric Brownian motion conditioned on endpoints:

```@example noise_types
# Bridge parameters
μ_bridge = 0.05
σ_bridge = 0.2
start_val = 100.0
end_val = 110.0

gbm_bridge = GeometricBrownianBridge(μ_bridge, σ_bridge, 0.0, start_val, 1.0, end_val)

prob_gbm_bridge = NoiseProblem(gbm_bridge, (0.0, 1.0))
sol_gbm_bridge = solve(prob_gbm_bridge; dt = 0.01)

println("GBM bridge starts at: $(sol_gbm_bridge.u[1])")
println("GBM bridge ends at: $(sol_gbm_bridge.u[end])")
```

## Process Comparisons

Let's compare different noise processes visually by looking at their statistical properties:

```@example noise_types
using Statistics

# Generate multiple realizations for statistical analysis
function analyze_process(noise_process, name, n_paths=1000)
    final_values = Float64[]
    
    for i in 1:n_paths
        # Reset the process for each simulation
        reset_process = typeof(noise_process)(noise_process.params...)
        prob = NoiseProblem(reset_process, (0.0, 1.0))
        sol = solve(prob; dt = 0.01)
        push!(final_values, sol.u[end])
    end
    
    println("$name Statistics (after t=1):")
    println("  Mean: $(round(mean(final_values), digits=4))")
    println("  Std:  $(round(std(final_values), digits=4))")
    println("  Min:  $(round(minimum(final_values), digits=4))")
    println("  Max:  $(round(maximum(final_values), digits=4))")
    println()
    
    return final_values
end

# Note: This analysis would require process constructors to be callable
# For demonstration, we'll just show the concept
println("Statistical analysis concept:")
println("- Wiener process: Mean ≈ 0, increasing variance with time")
println("- OU process: Mean reverts to long-term level")  
println("- GBM: Log-normal distribution of final values")
println("- Compound Poisson: Discrete jumps, possibly skewed distribution")
```

Each noise process has unique characteristics suitable for different modeling scenarios. Choose based on your specific application needs.