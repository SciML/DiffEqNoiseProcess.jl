module DiffEqNoiseProcessReverseDiffExt

using DiffEqNoiseProcess, DiffEqBase, Random
import ReverseDiff

@inline function DiffEqNoiseProcess.wiener_randn(
        rng::Random.AbstractRNG,
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
    return rand_vec .= ReverseDiff.track.(randn.((rng,), typeof.(DiffEqBase.value.(rand_vec))))
end
@inline function DiffEqNoiseProcess.wiener_randn!(
        rng::AbstractRNG,
        rand_vec::AbstractArray{
            <:ReverseDiff.TrackedReal,
        }
    )
    return rand_vec .= ReverseDiff.track.(randn.((rng,), typeof.(DiffEqBase.value.(rand_vec))))
end

end
