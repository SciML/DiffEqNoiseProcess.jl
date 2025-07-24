# DiffEqNoiseProcess.jl: Noise Processes for Stochastic Modeling

DiffEqNoiseProcess.jl provides a comprehensive suite of noise processes for stochastic differential equations and random differential equations. The `NoiseProcess` types are distributionally-exact, meaning they are generated directly according to their analytical distributions rather than as solutions to SDEs.

## Key Features

- **Mathematically rigorous**: All processes maintain distributional exactness
- **Comprehensive collection**: Wiener processes, Ornstein-Uhlenbeck, geometric Brownian motion, compound Poisson, and more
- **Flexible interface**: Works seamlessly with the DifferentialEquations.jl ecosystem
- **Memory efficient**: Includes virtual processes and wrappers for large-scale simulations
- **Bridging support**: Conditional processes for specific boundary conditions

## Installation

To install DiffEqNoiseProcess.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("DiffEqNoiseProcess")
```

## Quick Start

### Basic Wiener Process

Create and simulate a standard Brownian motion:

```@example index
using DiffEqNoiseProcess

# Create a Wiener process: WienerProcess(t0, W0, Z0)
W = WienerProcess(0.0, 0.0, 1.0)

# Simulate it over time interval [0, 1] 
prob = NoiseProblem(W, (0.0, 1.0))
sol = solve(prob; dt = 0.01)

println("Final Brownian motion value: $(sol.u[end])")
```

### Using Noise in SDE Problems

Noise processes integrate directly with SDE problems:

```@example index
using SciMLBase

# Define SDE: dX = μX dt + σX dW (geometric Brownian motion)
f(u, p, t) = 0.05 * u  # 5% drift
g(u, p, t) = 0.2 * u   # 20% volatility

# Use custom geometric Brownian motion noise
gbm_noise = GeometricBrownianMotionProcess(0.05, 0.2, 0.0, 100.0, 100.0)

# Create SDE problem (would normally solve this)
# prob = SDEProblem(f, g, 100.0, (0.0, 1.0), noise = gbm_noise)
```

### Ensemble Simulations

Generate multiple noise realizations:

```@example index
# Create ensemble problem for Monte Carlo simulation
enprob = EnsembleProblem(prob)
ensemble_sol = solve(enprob; dt = 0.01, trajectories = 100)

println("Generated $(length(ensemble_sol.u)) noise trajectories")
```

## Available Noise Processes

DiffEqNoiseProcess.jl includes a rich collection of noise processes:

### Classic Processes
- **WienerProcess**: Standard Brownian motion
- **RealWienerProcess**: Scalar Brownian motion  
- **CorrelatedWienerProcess**: Multi-dimensional correlated noise
- **OrnsteinUhlenbeckProcess**: Mean-reverting process
- **GeometricBrownianMotionProcess**: Financial modeling
- **CompoundPoissonProcess**: Jump processes

### Bridge Processes
- **BrownianBridge**: Brownian motion with fixed endpoints
- **GeometricBrownianBridge**: GBM with boundary conditions

### Advanced Features
- **NoiseFunction**: Custom noise from functions
- **NoiseGrid**: Noise from pre-computed data
- **NoiseWrapper**: Reuse previous simulations
- **VirtualBrownianTree**: Memory-efficient alternative
- **NoiseApproximation**: Colored noise via SDE solutions

## Documentation Structure

This documentation is organized as follows:

- **[Basic Usage Tutorial](tutorials/basic_usage.md)**: Get started with fundamental concepts
- **[Noise Process Types](tutorials/noise_processes.md)**: Comprehensive guide to all available processes
- **[Advanced Features](tutorials/advanced_features.md)**: Custom processes, wrappers, and optimization
- **[API Reference](api/noise_processes.md)**: Complete function and type documentation
- **[Interface API](api/interface.md)**: Step management and internal functions

## Mathematical Foundation

All noise processes in this package maintain mathematical rigor through:

1. **Distributional exactness**: Processes follow their theoretical distributions exactly
2. **Rejection Sampling with Memory (RSWM)**: Ensures correctness when adaptive stepping is used
3. **Proper bridging**: Conditional processes respect boundary conditions
4. **Interpolation consistency**: Values between grid points maintain distributional properties

## Contributing

  - Please refer to the
    [SciML ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://github.com/SciML/ColPrac/blob/master/README.md)
    for guidance on PRs, issues, and other matters relating to contributing to SciML.

  - See the [SciML Style Guide](https://github.com/SciML/SciMLStyle) for common coding practices and other style decisions.
  - There are a few community forums:
    
      + The #diffeq-bridged and #sciml-bridged channels in the
        [Julia Slack](https://julialang.org/slack/)
      + The #diffeq-bridged and #sciml-bridged channels in the
        [Julia Zulip](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
      + On the [Julia Discourse forums](https://discourse.julialang.org)
      + See also [SciML Community page](https://sciml.ai/community/)

## Reproducibility

```@raw html
<details><summary>The documentation of this SciML package was built using these direct dependencies,</summary>
```

```@example
using Pkg # hide
Pkg.status() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>and using this machine and Julia version.</summary>
```

```@example
using InteractiveUtils # hide
versioninfo() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>A more complete overview of all dependencies and their versions is also provided.</summary>
```

```@example
using Pkg # hide
Pkg.status(; mode = PKGMODE_MANIFEST) # hide
```

```@raw html
</details>
```

```@eval
using TOML
using Markdown
version = TOML.parse(read("../../Project.toml", String))["version"]
name = TOML.parse(read("../../Project.toml", String))["name"]
link_manifest = "https://github.com/SciML/" * name * ".jl/tree/gh-pages/v" * version *
                "/assets/Manifest.toml"
link_project = "https://github.com/SciML/" * name * ".jl/tree/gh-pages/v" * version *
               "/assets/Project.toml"
Markdown.parse("""You can also download the
[manifest]($link_manifest)
file and the
[project]($link_project)
file.
""")
```
