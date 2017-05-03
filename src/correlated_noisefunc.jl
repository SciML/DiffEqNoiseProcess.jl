function construct_correlated_noisefunc(Γ)
  γ = svdfact(Γ)
  A = γ[:U]*Diagonal(√γ[:S])
  b = Vector{eltype(Γ)}(size(Γ,1))
  dist = function (W,dt)
    if typeof(W.dW) <: AbstractArray
      return A*sqrt.(abs(dt))*wiener_randn(size(W.dW))
    else
      return A*sqrt.(abs(dt))*wiener_randn(typeof(W.dW))
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
  A = γ[:U]*Diagonal(√.(γ[:S]))
  b = Vector{eltype(Γ)}(size(Γ,1))
  dist = function (rand_vec,W,dt)
    wiener_randn!(b)
    b .*= sqrt.(abs(dt))
    A_mul_B!(rand_vec,A,b)
  end
  bridge = function (rand_vec,W,W0,Wh,q,h)
    error("Bridging distribution is unknown. Cannot use adapativity")
  end
  dist,bridge
end
CorrelatedWienerProcess!(Γ,t0,W0,Z0=nothing) = NoiseProcess(t0,W0,Z0,construct_correlated_noisefunc!(Γ)...,rswm=RSWM())
