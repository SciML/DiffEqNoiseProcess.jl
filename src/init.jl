wiener_randn!() = nothing

function __init__()
    @require CuArrays="3a865a2d-5b23-5a0f-bc46-62713ec82fae" begin
        wiener_randn!(rng::AbstractRNG,rand_vec::CuArrays.CuArray) = randn!(rand_vec)
    end
    
    @require CUDA="052768ef-5323-5732-b1bb-66c8b64840ba" begin
        wiener_randn!(rng::AbstractRNG,rand_vec::CUDA.CuArray) = randn!(rand_vec)
    end

    @require ReverseDiff="37e2e3b7-166d-5795-8a7a-e32c996b4267" begin
        @inline function DiffEqNoiseProcess.wiener_randn(rng::Random.AbstractRNG,proto::ReverseDiff.TrackedArray)
          ReverseDiff.track(convert.(eltype(proto.value),randn(rng,size(proto))))
        end
    end
end
