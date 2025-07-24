# Noise Processes API

## Wiener Processes

### Standard Wiener Process

```@docs
WienerProcess
WienerProcess!
SimpleWienerProcess
SimpleWienerProcess!
```

### Real-Valued Wiener Process

```@docs
RealWienerProcess
RealWienerProcess!
```

### Correlated Wiener Process

```@docs
CorrelatedWienerProcess
CorrelatedWienerProcess!
```

## Geometric Brownian Motion

```@docs
GeometricBrownianMotionProcess
GeometricBrownianMotionProcess!
```

## Ornstein-Uhlenbeck Process

```@docs
OrnsteinUhlenbeckProcess
OrnsteinUhlenbeckProcess!
```

## Jump Processes

```@docs
CompoundPoissonProcess
CompoundPoissonProcess!
```

## Bridge Processes

```@docs
BrownianBridge
BrownianBridge!
GeometricBrownianBridge
GeometricBrownianBridge!
OrnsteinUhlenbeckBridge
OrnsteinUhlenbeckBridge!
CompoundPoissonBridge
CompoundPoissonBridge!
```

## Advanced Noise Types

### Noise Wrapper

```@docs
NoiseWrapper
```

### Noise Functions

```@docs
NoiseFunction
NoiseTransport
```

### Noise from Data

```@docs
NoiseGrid
```

### Noise Approximation

```@docs
NoiseApproximation
```

### Memory-Efficient Alternatives

```@docs
VirtualBrownianTree
VirtualBrownianTree!
BoxWedgeTail
BoxWedgeTail!
```

## Preconditioned Crank-Nicolson

```@docs
pCN
pCN!
```