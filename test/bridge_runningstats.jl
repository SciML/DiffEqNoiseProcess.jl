# Shared streaming-statistics helpers for the bridge tests. Each @safetestset in
# bridge_test.jl runs in its own module, so this file is included into each one to
# provide the running mean/variance accumulator without retaining full solutions.
keep_u = (sol, ctx) -> ((u = sol.u,), false)

mutable struct RunningStats
    n::Int
    sum::Vector{Float64}
    sumsq::Vector{Float64}
    RunningStats() = new(0, Float64[], Float64[])
end
function push_stats!(rs::RunningStats, u)
    if rs.n == 0
        resize!(rs.sum, length(u))
        fill!(rs.sum, 0.0)
        resize!(rs.sumsq, length(u))
        fill!(rs.sumsq, 0.0)
    end
    rs.n += 1
    for i in eachindex(u)
        v = first(u[i])
        rs.sum[i] += v
        rs.sumsq[i] += abs2(v)
    end
    return rs
end
stats_mean(rs::RunningStats) = rs.sum ./ rs.n
function stats_std(rs::RunningStats)
    m = stats_mean(rs)
    return sqrt.(max.(rs.sumsq ./ rs.n .- abs2.(m), 0.0) .* (rs.n / (rs.n - 1)))
end
