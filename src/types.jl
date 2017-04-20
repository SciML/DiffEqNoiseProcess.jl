type NoiseProcess{T,T2,T3,ZType,F,F2,inplace,S1,S2,RSWM} <: AbstractNoiseProcess{inplace}
  dist::F
  bridge::F2
  t::Vector{T}
  W::Vector{T2}
  Z::ZType
  curt::T
  curW::T2
  curZ::T3
  dt::T
  dW::T2
  dZ::T3
  dWtilde::T2
  dZtilde::T3
  dWtmp::T2
  dZtmp::T3
  S₁::S1
  S₂::S2
  rswm::RSWM
  maxstacksize::Int
  maxstacksize2::Int
end
(W::NoiseProcess)(t) = interpolate!(W,t)
adaptive_alg(W::NoiseProcess) = adaptive_alg(W.rswm)
isinplace{T,T2,T3,ZType,F,F2,inplace,S1,S2,RSWM}(W::NoiseProcess{T,T2,T3,ZType,F,F2,inplace,S1,S2,RSWM}) = inplace

function NoiseProcess(t0,W0,Z0,dist,bridge;iip=DiffEqBase.isinplace(dist,3),
                       rswm = RSWM())
  S₁ = DataStructures.Stack{}(Tuple{typeof(t0),typeof(W0),typeof(Z0)})
  S₂ = ResettableStacks.ResettableStack{}(
                        Tuple{typeof(t0),typeof(W0),typeof(Z0)})
  if Z0==nothing
    Z=nothing
    curZ = nothing
    dZ = nothing
    dZtilde= nothing
    dZtmp = nothing
  else
    Z=[copy(Z0)]
    curZ = copy(Z0)
    dZ = copy(Z0)
    dZtilde= copy(Z0)
    dZtmp = copy(Z0)
  end
  NoiseProcess{typeof(t0),typeof(W0),typeof(dZ),typeof(Z),
                typeof(dist),typeof(bridge),
                iip,typeof(S₁),typeof(S₂),typeof(rswm)}(
                dist,bridge,[t0],[copy(W0)],Z,t0,
                copy(W0),curZ,t0,copy(W0),dZ,copy(W0),dZtilde,copy(W0),dZtmp,S₁,S₂,rswm,0,0)
end
