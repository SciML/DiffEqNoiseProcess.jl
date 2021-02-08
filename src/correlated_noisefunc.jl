function construct_correlated_noisefunc(Γ)
  γ = svd(Γ)
  A = γ.U*Diagonal(sqrt.(γ.S))
  dist = function (dW,W,dt,u,p,t,rng)
    if typeof(dW) <: AbstractArray
      return A*sqrt.(abs(dt))*wiener_randn(rng,dW)
    else
      return A*sqrt.(abs(dt))*wiener_randn(rng,typeof(dW))
    end
  end
  bridge = function (W,W0,Wh,q,h,u,p,t,rng)
    error("Bridging distribution is unknown. Cannot use adapativity")
  end
  dist,bridge
end
CorrelatedWienerProcess(Γ,t0,W0,Z0=nothing;rng = Xorshifts.Xoroshiro128Plus(rand(UInt64))) = NoiseProcess{false}(t0,W0,Z0,construct_correlated_noisefunc(Γ)...,rswm=RSWM(),rng=rng)

function construct_correlated_noisefunc!(Γ)
  γ = svd(Γ)
  A = γ.U*Diagonal(sqrt.(γ.S))
  b = Vector{eltype(Γ)}(undef,size(Γ,1))
  dist = function (rand_vec,W,dt,u,p,t,rng)
    wiener_randn!(rng,b)
    b .*= sqrt.(abs(dt))
    mul!(rand_vec,A,b)
  end
  bridge = function (rand_vec,W,W0,Wh,q,h,u,p,t,rng)
    error("Bridging distribution is unknown. Cannot use adapativity")
  end
  dist,bridge
end
CorrelatedWienerProcess!(Γ,t0,W0,Z0=nothing;rng = Xorshifts.Xoroshiro128Plus(rand(UInt64))) =
                         NoiseProcess{true}(t0,W0,Z0,construct_correlated_noisefunc!(Γ)...,
                         rswm=RSWM(),rng = rng)
