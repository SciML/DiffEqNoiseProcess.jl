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

function WHITE_NOISE_DIST(W,dt)
  if typeof(W.dW) <: AbstractArray
    return sqrt(abs(dt))*wiener_randn(size(W.dW))
  else
    return sqrt(abs(dt))*wiener_randn(typeof(W.dW))
  end
end
function WHITE_NOISE_BRIDGE(W,W0,Wh,q,h)
  if typeof(W.dW) <: AbstractArray
    return sqrt((1-q)*q*abs(h))*wiener_randn(size(W.dW))+q*(Wh-W0)+W0
  else
    return sqrt((1-q)*q*abs(h))*wiener_randn(typeof(W.dW))+q*(Wh-W0)+W0
  end
end
WienerProcess(t0,W0,Z0=nothing) = NoiseProcess(t0,W0,Z0,WHITE_NOISE_DIST,WHITE_NOISE_BRIDGE,rswm=RSWM())

function INPLACE_WHITE_NOISE_DIST(rand_vec,W,dt)
  wiener_randn!(rand_vec)
  rand_vec .*= sqrt(abs(dt))
end
function INPLACE_WHITE_NOISE_BRIDGE(rand_vec,W,W0,Wh,q,h)
  wiener_randn!(rand_vec)
  rand_vec .= sqrt((1.-q).*q.*abs(h)).*rand_vec.+q.*(Wh.-W0).+W0
end
WienerProcess!(t0,W0,Z0=nothing) = NoiseProcess(t0,W0,Z0,INPLACE_WHITE_NOISE_DIST,INPLACE_WHITE_NOISE_BRIDGE,rswm=RSWM())
