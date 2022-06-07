# Classic Noise Processes

This section describes the available `NoiseProcess` types. Note that all
keyword arguments are splatted into the `NoiseProcess` constructor, and thus
options like `reset` are available on the pre-built processes.

```@docs
WienerProcess
WienerProcess!
RealWienerProcess
RealWienerProcess!
OrnsteinUhlenbeckProcess
OrnsteinUhlenbeckProcess!
GeometricBrownianMotionProcess
GeometricBrownianMotionProcess!
CorrelatedWienerProcess
CorrelatedWienerProcess!
SimpleWienerProcess
SimpleWienerProcess!
```

## Bridges

```@docs
BrownianBridge
BrownianBridge!
GeometricBrownianBridge
GeometricBrownianBridge!
CompoundPoissonBridge
CompoundPoissonBridge!