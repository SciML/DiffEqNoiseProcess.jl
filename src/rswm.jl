"""
    RSWM(; discard_length = 1e-15, adaptivealg = :RSwM3)

Rejection Sampling with Memory (RSWM) algorithm configuration for noise processes.

RSWM ensures distributional exactness when adaptive time stepping is used with noise processes.
It maintains memory of rejected values to avoid biasing the noise distribution.

## Fields
- `discard_length`: Threshold for discarding stored values to save memory. Smaller values use more memory but are more accurate.
- `adaptivealg`: The adaptive algorithm variant to use (`:RSwM3` is recommended for Brownian motion)

## Algorithm Variants
- `:RSwM1`: Basic rejection sampling with single stack
- `:RSwM2`: Improved version with better memory management
- `:RSwM3`: Most advanced version with two-stack system (recommended for Brownian motion)
- `:RSwM0`: No memory storage variant for state-dependent noise processes

### About `:RSwM0` (No Memory)
For tau-leaping with state-dependent rates, the rate λ is approximated as constant over each
step. When a step is rejected and you use a Poisson bridge to pull back to a smaller dt,
storing the "future" portion doesn't make sense because you don't know the correct λ anyway -
it's always an approximation. Both storing and discarding have errors, but recalculating
fresh is more appropriate for this use case.

`:RSwM0` uses bridging/interpolation for rejected steps but discards future values instead
of storing them. This is the appropriate choice for:
- Compound Poisson processes with state-dependent rates
- Tau-leaping with post-leap adaptivity (Anderson's algorithm)

Reference: Anderson, D.F. "Incorporating postleap checks in tau-leaping"
J. Chem. Phys. 128, 054103 (2008); https://doi.org/10.1063/1.2819665

## Examples
```julia
# Conservative (high accuracy) for Brownian motion
rswm_accurate = RSWM(discard_length = 1e-12)

# Aggressive (lower memory usage) for Brownian motion
rswm_fast = RSWM(discard_length = 1e-6)

# No memory for state-dependent Poisson processes
rswm_nomem = RSWM(adaptivealg = :RSwM0)

# Use in noise process
W = WienerProcess(0.0, 0.0, 1.0; rswm = rswm_accurate)
```
"""
mutable struct RSWM{T}
    discard_length::T
    adaptivealg::Symbol
end

Base.@pure function RSWM(;
        discard_length = 1.0e-15,
        adaptivealg::Symbol = :RSwM3
    )
    RSWM{typeof(discard_length)}(discard_length, adaptivealg)
end

"""
    adaptive_alg(rswm::RSWM)

Get the adaptive algorithm type from an RSWM configuration.
"""
adaptive_alg(rswm::RSWM{T}) where {T} = rswm.adaptivealg
