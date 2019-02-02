isinplace(W::AbstractNoiseProcess{T,N,inplace}) where {T,N,inplace} = inplace

mutable struct NoiseProcess{T,N,Tt,T2,T3,ZType,F,F2,inplace,S1,S2,RSWM,RNGType} <: AbstractNoiseProcess{T,N,inplace}
  dist::F
  bridge::F2
  t::Vector{Tt}
  u::Vector{T2} # Aliased pointer to W for the AbstractVectorOfArray interface
  W::Vector{T2}
  Z::ZType
  curt::Tt
  curW::T2
  curZ::T3
  dt::Tt
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
  save_everystep::Bool
  iter::Int
  rng::RNGType
  reset::Bool
  reseed::Bool

  function NoiseProcess{iip}(t0,W0,Z0,dist,bridge;
                         rswm = RSWM(),save_everystep=true,
                         rng = Xorshifts.Xoroshiro128Plus(rand(UInt64)),
                         reset = true, reseed = true) where iip
    S₁ = DataStructures.Stack{Tuple{typeof(t0),typeof(W0),typeof(Z0)}}()
    S₂ = ResettableStacks.ResettableStack{iip}(
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
    W = [copy(W0)]
    N = length((size(W0)..., length(W)))
    new{eltype(eltype(W0)),N,typeof(t0),typeof(W0),typeof(dZ),typeof(Z),
                  typeof(dist),typeof(bridge),
                  iip,typeof(S₁),typeof(S₂),typeof(rswm),typeof(rng)}(
                  dist,bridge,[t0],W,W,Z,t0,
                  copy(W0),curZ,t0,copy(W0),dZ,copy(W0),dZtilde,copy(W0),dZtmp,
                  S₁,S₂,rswm,0,0,save_everystep,0,rng,reset,reseed)
  end

end
(W::NoiseProcess)(t) = interpolate!(W,t)
(W::NoiseProcess)(out1,out2,t) = interpolate!(out1,out2,W,t)
adaptive_alg(W::NoiseProcess) = adaptive_alg(W.rswm)

function NoiseProcess(t0,W0,Z0,dist,bridge;kwargs...)
  iip=DiffEqBase.isinplace(dist,4)
  NoiseProcess{iip}(t0,W0,Z0,dist,bridge;kwargs...)
end

mutable struct NoiseWrapper{T,N,Tt,T2,T3,T4,ZType,inplace} <: AbstractNoiseProcess{T,N,inplace}
  t::Vector{Tt}
  u::Vector{T2}
  W::Vector{T2}
  Z::ZType
  curt::Tt
  curW::T2
  curZ::T3
  dt::Tt
  dW::T2
  dZ::T3
  source::T4
  reset::Bool
end

function NoiseWrapper(source::AbstractNoiseProcess{T,N,inplace};
                      reset=true) where {T,N,inplace}
  if source.Z==nothing
    Z=nothing
    curZ = nothing
    dZ = nothing
  else
    Z=[copy(source.Z[1])]
    curZ = copy(source.Z[1])
    dZ = copy(source.Z[1])
  end
  W = [copy(source.W[1])]
  NoiseWrapper{T,N,typeof(source.t[1]),typeof(source.W[1]),typeof(dZ),typeof(source),typeof(Z),inplace}(
                [source.t[1]],W,W,Z,source.t[1],copy(source.W[1]),curZ,source.t[1],copy(source.W[1]),dZ,source,reset)
end

(W::NoiseWrapper)(t) = interpolate!(W,t)
(W::NoiseWrapper)(out1,out2,t) = interpolate!(out1,out2,W,t)
adaptive_alg(W::NoiseWrapper) = adaptive_alg(W.source)

mutable struct NoiseFunction{T,N,wType,zType,Tt,T2,T3,inplace} <: AbstractNoiseProcess{T,N,inplace}
  W::wType
  Z::zType
  curt::Tt
  curW::T2
  curZ::T3
  dt::Tt
  dW::T2
  dZ::T3
  reset::Bool

  function NoiseFunction{iip}(t0,W,Z=nothing;
                         noise_prototype=W(t0),reset=true) where iip
    curt = t0
    dt = t0
    curW = copy(noise_prototype)
    dW = copy(noise_prototype)
    if Z==nothing
      curZ = nothing
      dZ = nothing
    else
      curZ = copy(noise_prototype)
      dZ = copy(noise_prototype)
    end
    new{typeof(noise_prototype),ndims(noise_prototype),typeof(W),typeof(Z),
                  typeof(curt),typeof(curW),typeof(curZ),iip}(W,Z,curt,curW,curZ,
                  dt,dW,dZ,reset)
  end

end

function (W::NoiseFunction)(t)
  if W.Z != nothing
    if isinplace(W)
      out2 = similar(W.dZ)
      W.Z(out2,t)
    else
      out2 = W.Z(t)
    end
  else
    out2 = nothing
  end
  if isinplace(W)
    out1 = similar(W.dW)
    W.W(out1,t)
  else
    out1 = W.W(t)
  end
  out1,out2
end
function (W::NoiseFunction)(out1,out2,t)
  W.W(out1,t)
  W.Z != nothing && W.Z(out2,t)
end

function NoiseFunction(t0,W,Z=nothing;kwargs...)
  iip=DiffEqBase.isinplace(W,2)
  NoiseFunction{iip}(t0,W,Z;kwargs...)
end

mutable struct NoiseGrid{T,N,Tt,T2,T3,ZType,inplace} <: AbstractNoiseProcess{T,N,inplace}
  t::Vector{Tt}
  u::Vector{T2}
  W::Vector{T2}
  Z::ZType
  curt::Tt
  curW::T2
  curZ::T3
  dt::Tt
  dW::T2
  dZ::T3
  step_setup::Bool
  reset::Bool
end

function NoiseGrid(t,W,Z=nothing;reset=true)
  val = W[1]
  curt = t[1]
  dt = t[1]
  curW = copy(val)
  dW = copy(val)
  if Z==nothing
    curZ = nothing
    dZ = nothing
  else
    curZ = copy(Z[1])
    dZ = copy(Z[1])
  end
  typeof(val) <: AbstractArray ? iip = true : iip = false
  NoiseGrid{typeof(val),ndims(val),typeof(dt),typeof(dW),typeof(dZ),typeof(Z),iip}(
            t,W,W,Z,curt,curW,curZ,dt,dW,dZ,true,reset)
end

(W::NoiseGrid)(t) = interpolate!(W,t)
(W::NoiseGrid)(out1,out2,t) = interpolate!(out1,out2,W,t)

mutable struct NoiseApproximation{T,N,Tt,T2,T3,S1,S2,ZType,inplace} <: AbstractNoiseProcess{T,N,inplace}
  t::Vector{Tt}
  u::Vector{T2}
  W::Vector{T2}
  Z::ZType
  curt::Tt
  curW::T2
  curZ::T3
  dt::Tt
  dW::T2
  dZ::T3
  source1::S1
  source2::S2
  reset::Bool
end

function NoiseApproximation(source1::DEIntegrator,source2::Union{DEIntegrator,Nothing}=nothing;
                   reset=true)
  _source1 = deepcopy(source1)
  _source2 = deepcopy(source2)
  if _source2==nothing
    Z=nothing
    curZ = nothing
    dZ = nothing
  else
    Z=_source2.sol.u
    curZ = copy(_source2.u)
    dZ = copy(_source2.u)
    _source2.opts.advance_to_tstop = true
  end
  val = copy(_source1.u)
  t = _source1.sol.t
  W = _source1.sol.u
  curW = copy(_source1.u)
  dW = copy(_source1.u)
  dt = _source1.dt
  curt = _source1.t
  _source1.opts.advance_to_tstop = true
  NoiseApproximation{typeof(val),ndims(val),typeof(curt),typeof(curW),typeof(curZ),
                     typeof(_source1),typeof(_source2),typeof(Z),
                     isinplace(_source1.sol.prob)}(
                     t,W,W,Z,curt,curW,curZ,dt,dW,dZ,_source1,_source2,reset)
end

(W::NoiseApproximation)(t) = interpolate!(W,t)
(W::NoiseApproximation)(out1,out2,t) = interpolate!(out1,out2,W,t)
