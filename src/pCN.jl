# Preconditioned Crank–Nicolson algorithm tools

function generate_innovation(W0::Number, t, rng)
    dt = diff(t)
    Wnew = cumsum([zero(W0);
                   [sqrt(dt[i]) .* wiener_randn(rng, typeof(W0))
                    for i in 1:(length(t) - 1)]])
end

function generate_innovation(W0, t, rng)
    dt = diff(t)
    Wnew = cumsum([[zero(W0)];
                   [sqrt(dt[i]) .* wiener_randn(rng, W0)
                    for i in 1:(length(t) - 1)]])
end

@doc doc"""
    pCN!(noise::AbstractNoiseProcess, ρ; reset=true,reverse=false,indx=nothing)

Create a new, but correlated noise process from `noise` and additional entropy with correlation ρ.
This update defines an autoregressive process in the space of Wiener (or noise process) trajectories, which can be used as proposal distribution in Metropolis-Hastings algorithms (often called the "preconditioned Crank–Nicolson scheme".)

External links
 * [Preconditioned Crank–Nicolson algorithm on Wikipedia](https://en.wikipedia.org/wiki/Preconditioned_Crank–Nicolson_algorithm)
"""
function pCN!(source::AbstractNoiseProcess{T, N, Vector{T2}, inplace}, ρ;
        reset = true, reverse = false, indx = nothing) where {T, N, T2, inplace}

    # generate new Wiener process similar to the one in source
    Wnew = generate_innovation(source.W[1], source.t, source.rng)

    source.W = ρ * source.W + sqrt(one(ρ) - ρ^2) * Wnew
    source.u = ρ * source.u + sqrt(one(ρ) - ρ^2) * Wnew
    NoiseWrapper(source, reset = reset, reverse = reverse, indx = indx)
end

"""
    pCN(noise::AbstractNoiseProcess, ρ; reset=true,reverse=false,indx=nothing)

Create a new, but correlated noise process from `noise` and additional entropy with correlation ρ.
"""
function pCN(source::AbstractNoiseProcess{T, N, Vector{T2}, inplace}, ρ;
        reset = true, reverse = false, indx = nothing) where {T, N, T2, inplace}
    source′ = deepcopy(source)

    # generate new Wiener process similar to the one in source
    Wnew = generate_innovation(source′.W[1], source′.t, source′.rng)

    source′.W = ρ * source′.W + sqrt(one(ρ) - ρ^2) * Wnew
    source′.u = ρ * source′.u + sqrt(one(ρ) - ρ^2) * Wnew
    NoiseWrapper(source′, reset = reset, reverse = reverse, indx = indx)
end

"""
    pCN(noise::NoiseGrid, ρ; reset=true, rng = Xorshifts.Xoroshiro128Plus(rand(UInt64)))

Create a new, but correlated noise process from `noise` and additional entropy with correlation ρ.
This update defines an autoregressive process in the space of Wiener (or noise process) trajectories, which can be used as proposal distribution in Metropolis-Hastings algorithms (often called the "preconditioned Crank–Nicolson scheme".)

External links

  - [Preconditioned Crank–Nicolson algorithm on Wikipedia](https://en.wikipedia.org/wiki/Preconditioned_Crank%E2%80%93Nicolson_algorithm)
"""
function pCN(source::NoiseGrid, ρ; reset = true,
        rng = Xorshifts.Xoroshiro128Plus(rand(UInt64)))

    # generate new Wiener process similar to the one in source
    t = source.t
    Wnew = generate_innovation(source.W[1], source.t, rng)

    W = ρ * source.W + sqrt(one(ρ) - ρ^2) * Wnew
    NoiseGrid(t, W, source.Z; reset = reset)
end
