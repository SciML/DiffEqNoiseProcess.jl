module DiffEqNoiseProcessReverseDiffExt

using DiffEqNoiseProcess: DiffEqNoiseProcess
using Random: AbstractRNG
using SciMLBase: value
import ReverseDiff

# TrackedArray storage cannot be mutated (no setindex!/broadcast-destination
# since ReverseDiff 1.18.6), while in-place noise processes require mutable
# buffers. collect gives a plain Array{TrackedReal} whose elements keep their
# tape origin, matching the Array{TrackedReal} convention used for u0.
# Out-of-place processes keep TrackedArray storage: they rebind rather than
# mutate, and wiener_randn on a TrackedArray prototype produces the draws.
@inline function DiffEqNoiseProcess.NoiseProcess{true}(
        t0, W0::ReverseDiff.TrackedArray, Z0, dist, bridge; kwargs...
    )
    return DiffEqNoiseProcess.NoiseProcess{true}(
        t0, collect(W0),
        Z0 isa ReverseDiff.TrackedArray ? collect(Z0) : Z0,
        dist, bridge; kwargs...
    )
end

@inline function DiffEqNoiseProcess.wiener_randn(
        rng::AbstractRNG,
        proto::ReverseDiff.TrackedArray
    )
    return ReverseDiff.track(convert.(eltype(proto.value), randn(rng, size(proto))))
end
@inline function DiffEqNoiseProcess.wiener_randn!(
        rng::AbstractRNG,
        rand_vec::Array{
            <:ReverseDiff.TrackedReal,
        }
    )
    return rand_vec .= ReverseDiff.track.(randn.((rng,), typeof.(value.(rand_vec))))
end
@inline function DiffEqNoiseProcess.wiener_randn!(
        rng::AbstractRNG,
        rand_vec::AbstractArray{
            <:ReverseDiff.TrackedReal,
        }
    )
    return rand_vec .= ReverseDiff.track.(randn.((rng,), typeof.(value.(rand_vec))))
end

end
