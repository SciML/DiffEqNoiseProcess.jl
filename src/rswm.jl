"""
    RSWM(; discard_length = 1e-15, adaptivealg = :RSwM3)

Rejection Sampling with Memory (RSWM) algorithm configuration for noise processes.

RSWM ensures distributional exactness when adaptive time stepping is used with noise processes.
It maintains memory of rejected values to avoid biasing the noise distribution.

## Fields
- `discard_length`: Threshold for discarding stored values to save memory. Smaller values use more memory but are more accurate.
- `adaptivealg`: The adaptive algorithm variant to use (`:RSwM3` is recommended)

## Examples
```julia
# Conservative (high accuracy)
rswm_accurate = RSWM(discard_length = 1e-12)

# Aggressive (lower memory usage)  
rswm_fast = RSWM(discard_length = 1e-6)

# Use in noise process
W = WienerProcess(0.0, 0.0, 1.0; rswm = rswm_accurate)
```
"""
mutable struct RSWM{T}
    discard_length::T
    adaptivealg::Symbol
end

Base.@pure function RSWM(;
        discard_length = 1e-15,
        adaptivealg::Symbol = :RSwM3)
    RSWM{typeof(discard_length)}(discard_length, adaptivealg)
end

"""
    adaptive_alg(rswm::RSWM)

Get the adaptive algorithm type from an RSWM configuration.
"""
adaptive_alg(rswm::RSWM{T}) where {T} = rswm.adaptivealg
