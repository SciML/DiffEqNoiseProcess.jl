function construct_correlated_noisefunc(Γ)
  γ = svdfact(Γ)
  A = γ[:U]*Diagonal(√γ[:S])
  b = Vector{eltype(Γ)}(size(Γ,1))
  dist = function (W,dt)
    if typeof(W.dW) <: AbstractArray
      return A*sqrt(abs(dt))*wiener_randn(size(W.dW))
    else
      return A*sqrt(abs(dt))*wiener_randn(typeof(W.dW))
    end
  end
  bridge = function (W,W0,Wh,q,h)
    error("Bridging distribution is unknown. Cannot use adapativity")
  end
  dist,bridge
end
CorrelatedWienerProcess(Γ,t0,W0,Z0=nothing) = NoiseProcess(t0,W0,Z0,construct_correlated_noisefunc(Γ)...,rswm=RSWM())

function construct_correlated_noisefunc!(Γ)
  γ = svdfact(Γ)
  A = γ[:U]*Diagonal(√γ[:S])
  b = Vector{eltype(Γ)}(size(Γ,1))
  dist = function (rand_vec,W,dt)
    wiener_randn!(b)
    b .*= sqrt(abs(dt))
    A_mul_B!(rand_vec,A,b)
  end
  bridge = function (rand_vec,W,W0,Wh,q,h)
    error("Bridging distribution is unknown. Cannot use adapativity")
  end
  dist,bridge
end
CorrelatedWienerProcess!(Γ,t0,W0,Z0=nothing) = NoiseProcess(t0,W0,Z0,construct_correlated_noisefunc!(Γ)...,rswm=RSWM())

"""
construct_correlated_noisefunc(Γ::AbstractArray)

Takes in a constant Covariance matrix Γ and spits out the noisefunc.
"""
function construct_correlated_noisefunc(Γ::AbstractArray)
  γ = svdfact(Γ)
  A = γ[:U]*Diagonal(√γ[:S])
  b = Vector{eltype(Γ)}(size(Γ,1))
  noise_func! = function (a,integrator)
    randn!(b)
    A_mul_B!(a,A,b)
  end
  NoiseProcess{:White,true,typeof(noise_func!)}(noise_func!)
end
