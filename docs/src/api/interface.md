# Noise Process Interface

This page documents the interface functions for working with noise processes.

## Step Management

```@docs
accept_step!
reject_step!
calculate_step!
setup_next_step!
save_noise!
```

## Noise Process Types

### Abstract Types

```@docs
AbstractNoiseProcess
```

### Core Types

```@docs
NoiseProcess
SimpleNoiseProcess
```

## Configuration

### Rejection Sampling with Memory (RSWM)

```@docs
RSWM
adaptive_alg
```

## Internal Functions

These functions are used internally by the noise process implementations:

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
DiffEqNoiseProcess.CompoundPoissonProcess
DiffEqNoiseProcess.CompoundPoissonProcess!
DiffEqNoiseProcess.cpp_bridge
DiffEqNoiseProcess.cpp_bridge!
```