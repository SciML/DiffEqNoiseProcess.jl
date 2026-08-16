# Preconditioned Crank–Nicolson algorithm tools

function generate_innovation(W0::Number, t, rng)
    dt = diff(t)
    return Wnew = cumsum(
        [
            zero(W0);
            [
                sqrt(dt[i]) .* wiener_randn(rng, typeof(W0))
                    for i in 1:(length(t) - 1)
            ]
        ]
    )
end

function generate_innovation(W0, t, rng)
    dt = diff(t)
    return Wnew = cumsum(
        [
            [zero(W0)];
            [
                sqrt(dt[i]) .* wiener_randn(rng, W0)
                    for i in 1:(length(t) - 1)
            ]
        ]
    )
end

"""
    pCN!(noise::AbstractNoiseProcess, ρ; reset=true,reverse=false,indx=nothing)

Create a correlated noise proposal in-place using the preconditioned
Crank–Nicolson update. The source path is replaced by
`ρ * source + sqrt(1 - ρ^2) * innovation`.

# Arguments
- `noise`: Source `AbstractNoiseProcess` whose saved path is updated.
- `ρ`: Correlation parameter, normally in `[-1, 1]`.

# Keywords
- `reset`: Whether the returned wrapper resets before solving.
- `reverse`: Whether the returned wrapper traverses the source backwards.
- `indx`: Initial source index for the wrapper; defaults to the start (or end
  when `reverse = true`).

# Returns
A `NoiseWrapper` around the updated source. The source itself is mutated.

# Example
```julia
W = WienerProcess(0.0, 0.0)
proposal = pCN!(W, 0.9)
```

External links
 * [Preconditioned Crank–Nicolson algorithm on Wikipedia](https://en.wikipedia.org/wiki/Preconditioned_Crank–Nicolson_algorithm)
"""
function pCN!(
        source::AbstractNoiseProcess{T, N, Vector{T2}, inplace}, ρ;
        reset = true, reverse = false, indx = nothing
    ) where {T, N, T2, inplace}

    # generate new Wiener process similar to the one in source
    Wnew = generate_innovation(source.W[1], source.t, source.rng)

    source.W = ρ * source.W + sqrt(one(ρ) - ρ^2) * Wnew
    source.u = ρ * source.u + sqrt(one(ρ) - ρ^2) * Wnew
    return NoiseWrapper(source, reset = reset, reverse = reverse, indx = indx)
end

"""
    pCN(noise::AbstractNoiseProcess, ρ; reset=true,reverse=false,indx=nothing)

Create a correlated noise proposal from a copy of `noise` using the
preconditioned Crank–Nicolson update. The input process is not mutated.

# Arguments
- `noise`: Source `AbstractNoiseProcess` with a saved path.
- `ρ`: Correlation parameter, normally in `[-1, 1]`.

# Keywords
- `reset`: Whether the returned wrapper resets before solving.
- `reverse`: Whether the returned wrapper traverses the source backwards.
- `indx`: Initial source index for the wrapper.

# Returns
A `NoiseWrapper` containing the correlated proposal.

# Example
```julia
W = WienerProcess(0.0, 0.0)
proposal = pCN(W, 0.9)
```
"""
function pCN(
        source::AbstractNoiseProcess{T, N, Vector{T2}, inplace}, ρ;
        reset = true, reverse = false, indx = nothing
    ) where {T, N, T2, inplace}
    source′ = copy(source)

    # generate new Wiener process similar to the one in source
    Wnew = generate_innovation(source′.W[1], source′.t, source′.rng)

    source′.W = ρ * source′.W + sqrt(one(ρ) - ρ^2) * Wnew
    source′.u = ρ * source′.u + sqrt(one(ρ) - ρ^2) * Wnew
    return NoiseWrapper(source′, reset = reset, reverse = reverse, indx = indx)
end

"""
    pCN(noise::NoiseGrid, ρ; reset=true, rng = Random.default_rng())

Create a correlated proposal from a `NoiseGrid`. The source grid is not
mutated; the returned grid has the same time points and auxiliary process.

# Arguments
- `noise`: Source `NoiseGrid`.
- `ρ`: Correlation parameter, normally in `[-1, 1]`.

# Keywords
- `reset`: Reset flag stored on the returned grid.
- `rng`: Random number generator used for the innovation.

# Returns
A new `NoiseGrid` containing the correlated proposal.

# Example
```julia
grid = NoiseGrid(0.0:0.1:1.0, zeros(11))
proposal = pCN(grid, 0.9)
```

External links

  - [Preconditioned Crank–Nicolson algorithm on Wikipedia](https://en.wikipedia.org/wiki/Preconditioned_Crank%E2%80%93Nicolson_algorithm)
"""
function pCN(
        source::NoiseGrid, ρ; reset = true,
        rng = Random.default_rng()
    )

    # generate new Wiener process similar to the one in source
    t = source.t
    Wnew = generate_innovation(source.W[1], source.t, rng)

    W = ρ * source.W + sqrt(one(ρ) - ρ^2) * Wnew
    return NoiseGrid(t, W, source.Z; reset = reset)
end
