const one_over_sqrt2 = 1/sqrt(2)
@inline wiener_randn(rng::AbstractRNG,::Type{T}) where T = randn(rng,T)
@inline function wiener_randn(rng::AbstractRNG,proto::Array{T}) where T
  randn(rng,size(proto))
end
@inline function wiener_randn(rng::AbstractRNG,proto)
  convert(typeof(proto),randn(rng,size(proto)))
end
@inline wiener_randn!(rng::AbstractRNG,rand_vec::Array) = randn!(rng,rand_vec)

# TODO: This needs an overload for GPUs
@inline wiener_randn!(rng::AbstractRNG,rand_vec) = rand_vec .= Base.Broadcast.Broadcasted(randn,())
@inline wiener_randn(y::AbstractRNG,::Type{Complex{T}}) where T = convert(T,one_over_sqrt2)*(randn(y,T)+im*randn(y,T))

@inline function wiener_randn!(y::AbstractRNG,x::AbstractArray{<:Complex{T}}) where T<:Number
  # Remove loop
  @inbounds for i in eachindex(x)
    x[i] = convert(T,one_over_sqrt2)*(randn(y,T)+im*randn(y,T))
  end
end

@inline function WHITE_NOISE_DIST(W,dt,rng)
  if typeof(W.dW) <: AbstractArray && !(typeof(W.dW) <: SArray)
    return @fastmath sqrt(abs(dt))*wiener_randn(rng,W.dW)
  else
    return @fastmath sqrt(abs(dt))*wiener_randn(rng,typeof(W.dW))
  end
end

function WHITE_NOISE_BRIDGE(W,W0,Wh,q,h,rng)
  if typeof(W.dW) <: AbstractArray
    return @fastmath sqrt((1-q)*q*abs(h))*wiener_randn(rng,W.dW)+q*Wh
  else
    return @fastmath sqrt((1-q)*q*abs(h))*wiener_randn(rng,typeof(W.dW))+q*Wh
  end

end
WienerProcess(t0,W0,Z0=nothing;kwargs...) = NoiseProcess{false}(t0,W0,Z0,WHITE_NOISE_DIST,WHITE_NOISE_BRIDGE;kwargs...)

function INPLACE_WHITE_NOISE_DIST(rand_vec,W,dt,rng)
  wiener_randn!(rng,rand_vec)
  sqrtabsdt = @fastmath sqrt(abs(dt))
  @.. rand_vec *= sqrtabsdt
end
function INPLACE_WHITE_NOISE_BRIDGE(rand_vec,W,W0,Wh,q,h,rng)
  wiener_randn!(rng,rand_vec)
  #rand_vec .= sqrt((1.-q).*q.*abs(h)).*rand_vec.+q.*Wh
  sqrtcoeff = @fastmath sqrt((1-q)*q*abs(h))
  @.. rand_vec = sqrtcoeff*rand_vec+q*Wh
end
WienerProcess!(t0,W0,Z0=nothing;kwargs...) = NoiseProcess{true}(t0,W0,Z0,INPLACE_WHITE_NOISE_DIST,INPLACE_WHITE_NOISE_BRIDGE;kwargs...)



#### Real Valued Wiener Process. Ignores complex and the like
function REAL_WHITE_NOISE_DIST(W,dt,rng)
  if typeof(W.dW) <: AbstractArray
    return @fastmath sqrt(abs(dt))*randn(rng,size(W.dW))
  else
    return @fastmath sqrt(abs(dt))*randn(rng)
  end
end
function REAL_WHITE_NOISE_BRIDGE(W,W0,Wh,q,h,rng)
  if typeof(W.dW) <: AbstractArray
    return @fastmath sqrt((1-q)*q*abs(h))*randn(rng,size(W.dW))+q*Wh
  else
    return @fastmath sqrt((1-q)*q*abs(h))*randn(rng)+q*Wh
  end
end
RealWienerProcess(t0,W0,Z0=nothing;kwargs...) = NoiseProcess{false}(t0,W0,Z0,REAL_WHITE_NOISE_DIST,REAL_WHITE_NOISE_BRIDGE;kwargs...)

function REAL_INPLACE_WHITE_NOISE_DIST(rand_vec,W,dt,rng)
  sqabsdt = @fastmath sqrt(abs(dt))
  wiener_randn!(rng,rand_vec)
  sqabsdt = sqrt(abs(dt))
  @.. rand_vec *= sqabsdt
end
function REAL_INPLACE_WHITE_NOISE_BRIDGE(rand_vec,W,W0,Wh,q,h,rng)
  wiener_randn!(rng,rand_vec!)
  @.. rand_vec = @fastmath sqrt((1-q)*q*abs(h))*rand_vec+q*Wh
end
RealWienerProcess!(t0,W0,Z0=nothing;kwargs...) = NoiseProcess{true}(t0,W0,Z0,REAL_INPLACE_WHITE_NOISE_DIST,REAL_INPLACE_WHITE_NOISE_BRIDGE;kwargs...)
