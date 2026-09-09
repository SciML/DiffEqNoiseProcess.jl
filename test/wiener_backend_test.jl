using DiffEqNoiseProcess, Random, Test

struct BackendArray{T, N} <: AbstractArray{T, N}
    data::Array{T, N}
end

Base.size(array::BackendArray) = size(array.data)
Base.getindex(array::BackendArray, index::Int) = getindex(array.data, index)
Base.IndexStyle(::Type{<:BackendArray}) = IndexLinear()
const backend_randn_called = Ref(false)
function Base.similar(::BackendArray, ::Type{T}, dimensions::Dims) where {T}
    return BackendArray(Array{T}(undef, dimensions))
end
function Random.randn!(rng::AbstractRNG, array::BackendArray)
    backend_randn_called[] = true
    randn!(rng, array.data)
    return array
end

@testset "Wiener samples preserve the prototype backend" begin
    seed = UInt64(0x5eed)
    prototype = BackendArray(zeros(Float32, 2, 3))
    sample = DiffEqNoiseProcess.wiener_randn(Xoshiro(seed), prototype)

    @test backend_randn_called[]
    @test sample isa BackendArray
    @test sample.data == randn(Xoshiro(seed), Float32, size(prototype))
end
