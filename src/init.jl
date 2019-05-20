wiener_randn!() = nothing

function __init__()
    @require CuArrays="3a865a2d-5b23-5a0f-bc46-62713ec82fae" begin
        wiener_randn!(rng::AbstractRNG,rand_vec::CuArrays.CuArray) = randn!(rand_vec)
    end
end
