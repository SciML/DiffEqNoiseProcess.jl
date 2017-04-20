one_over_sqrt2 = 1/sqrt(2)
@inline wiener_randn() = randn()
@inline wiener_randn(x...) = randn(x...)
@inline wiener_randn!(x...) = randn!(x...)
@inline wiener_randn{T<:Number}(::Type{Complex{T}}) = one_over_sqrt2*(randn(T)+im*randn(T))
@inline wiener_randn{T<:Number}(::Type{Complex{T}},x...) = one_over_sqrt2*(randn(T,x...)+im*randn(T,x...))
@inline wiener_randn{T<:Number}(y::AbstractRNG,::Type{Complex{T}},x...) = one_over_sqrt2*(randn(y,T,x...)+im*randn(y,T,x...))
@inline wiener_randn{T<:Number}(y::AbstractRNG,::Type{Complex{T}}) = one_over_sqrt2*(randn(y,T)+im*randn(y,T))
@inline function wiener_randn!{T<:Number}(y::AbstractRNG,x::AbstractArray{Complex{T}})
  for i in eachindex(x)
    x[i] = one_over_sqrt2*(randn(y,T)+im*randn(y,T))
  end
end
@inline function wiener_randn!{T<:Number}(x::AbstractArray{Complex{T}})
  for i in eachindex(x)
    x[i] = one_over_sqrt2*(randn(T)+im*randn(T))
  end
end

WHITE_NOISE_DIST  = function (W,dt)
  if typeof(W.dW) <: AbstractArray
    return sqrt(dt)*wiener_randn(size(W.dW))
  else
    return sqrt(dt)*wiener_randn(typeof(W.dW))
  end
end
WHITE_NOISE_BRIDGE= (W,W0,Wh,q,h) -> sqrt((1-q)*q*h)*wiener_randn(typeof(W.dW))+q*(Wh-W0)+W0

WienerProcess(t0,W0) = NoiseProcess(t0,W0,WHITE_NOISE_DIST,WHITE_NOISE_BRIDGE,rswm=RSWM())

INPLACE_WHITE_NOISE_DIST  = function (rand_vec,W,dt)
  wiener_randn!(rand_vec)
  rand_vec .*= sqrt(dt)
end
INPLACE_WHITE_NOISE_BRIDGE = function (rand_vec,W,W0,Wh,q,h)
  wiener_randn!(rand_vec)
  rand_vec .= sqrt((1.-q).*q.*h).*rand_vec.+q.*(Wh.-W0).+W0
end
WienerProcess!(t0,W0) = NoiseProcess(t0,W0,INPLACE_WHITE_NOISE_DIST,INPLACE_WHITE_NOISE_BRIDGE,rswm=RSWM())
