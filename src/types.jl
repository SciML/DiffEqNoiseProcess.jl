isinplace(W::AbstractNoiseProcess{T,N,S,inplace}) where {T,N,S,inplace} = inplace

mutable struct NoiseProcess{T,N,Tt,T2,T3,ZType,F,F2,inplace,S1,S2,RSWM,C,RNGType} <: AbstractNoiseProcess{T,N,Vector{T2},inplace}
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
  continuous::Bool
  cache::C

  function NoiseProcess{iip}(t0,W0,Z0,dist,bridge;
                         rswm = RSWM(),save_everystep=true,
                         rng = Xorshifts.Xoroshiro128Plus(rand(UInt64)),
                         reset = true, reseed = true, continuous = true,
                         cache = nothing) where iip
    S₁ = ResettableStacks.ResettableStack{iip}(
                          Tuple{typeof(t0),typeof(W0),typeof(Z0)})
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
                  iip,typeof(S₁),typeof(S₂),typeof(rswm),typeof(cache),typeof(rng)}(
                  dist,bridge,[t0],W,W,Z,t0,
                  copy(W0),curZ,t0,copy(W0),dZ,copy(W0),dZtilde,copy(W0),dZtmp,
                  S₁,S₂,rswm,0,0,save_everystep,0,rng,reset,reseed,continuous,cache)
  end

end
(W::NoiseProcess)(t) = interpolate!(W,nothing,nothing,t)
(W::NoiseProcess)(u,p,t) = interpolate!(W,u,p,t)
(W::NoiseProcess)(out1,out2,u,p,t) = interpolate!(out1,out2,W,u,p,t)
adaptive_alg(W::NoiseProcess) = adaptive_alg(W.rswm)

function NoiseProcess(t0,W0,Z0,dist,bridge;kwargs...)
  iip=DiffEqBase.isinplace(dist,7)
  NoiseProcess{iip}(t0,W0,Z0,dist,bridge;kwargs...)
end

mutable struct SimpleNoiseProcess{T,N,Tt,T2,T3,ZType,F,F2,inplace,RNGType} <: AbstractNoiseProcess{T,N,Vector{T2},inplace}
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
  save_everystep::Bool
  iter::Int
  rng::RNGType
  reset::Bool
  reseed::Bool

  function SimpleNoiseProcess{iip}(t0,W0,Z0,dist,bridge;
                         save_everystep=true,
                         rng = Xorshifts.Xoroshiro128Plus(rand(UInt64)),
                         reset = true, reseed = true) where iip
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
                  typeof(dist),typeof(bridge),iip,typeof(rng)}(
                  dist,bridge,[t0],W,W,Z,t0,
                  copy(W0),curZ,t0,copy(W0),dZ,copy(W0),dZtilde,copy(W0),dZtmp,
                  save_everystep,0,rng,reset,reseed)
  end

end
(W::SimpleNoiseProcess)(t) = interpolate!(W,nothing,nothing,t)
(W::SimpleNoiseProcess)(u,p,t) = interpolate!(W,u,p,t)
(W::SimpleNoiseProcess)(out1,out2,u,p,t) = interpolate!(out1,out2,W,u,p,t)

function SimpleNoiseProcess(t0,W0,Z0,dist,bridge;kwargs...)
  iip=DiffEqBase.isinplace(dist,7)
  SimpleNoiseProcess{iip}(t0,W0,Z0,dist,bridge;kwargs...)
end

mutable struct NoiseWrapper{T,N,Tt,T2,T3,T4,ZType,inplace} <: AbstractNoiseProcess{T,N,Vector{T2},inplace}
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
  reverse::Bool
end

function NoiseWrapper(source::AbstractNoiseProcess{T,N,Vector{T2},inplace};
                      reset=true,reverse=false,indx=nothing) where {T,N,T2,inplace}

  if indx === nothing
    if reverse
      indx = length(source.t)
    else
      indx = 1
    end
  end

  if source.Z==nothing
    Z=nothing
    curZ = nothing
    dZ = nothing
  else
    Z=[copy(source.Z[indx])]
    curZ = copy(source.Z[indx])
    dZ = copy(source.Z[indx])
  end
  W = [copy(source.W[indx])]

  NoiseWrapper{T,N,typeof(source.t[1]),typeof(source.W[1]),typeof(dZ),typeof(source),typeof(Z),inplace}(
                [source.t[indx]],W,W,Z,source.t[indx],copy(source.W[indx]),curZ,source.t[indx],copy(source.W[indx]),dZ,source,reset,reverse)
end

(W::NoiseWrapper)(t) = interpolate!(W,nothing,nothing,t)
(W::NoiseWrapper)(u,p,t) = interpolate!(W,u,p,t)
(W::NoiseWrapper)(out1,out2,u,p,t) = interpolate!(out1,out2,W,u,p,t,reverse=W.reverse)
adaptive_alg(W::NoiseWrapper) = adaptive_alg(W.source)

mutable struct NoiseFunction{T,N,wType,zType,Tt,T2,T3,inplace} <: AbstractNoiseProcess{T,N,nothing,inplace}
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
                         noise_prototype=W(nothing,nothing,t0),reset=true) where iip
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

(W::NoiseFunction)(t) = W(nothing,nothing,t)
function (W::NoiseFunction)(u,p,t)
  if W.Z != nothing
    if isinplace(W)
      out2 = similar(W.dZ)
      W.Z(out2,u,p,t)
    else
      out2 = W.Z(u,p,t)
    end
  else
    out2 = nothing
  end
  if isinplace(W)
    out1 = similar(W.dW)
    W.W(out1,u,p,t)
  else
    out1 = W.W(u,p,t)
  end
  out1,out2
end
function (W::NoiseFunction)(out1,out2,u,p,t)
  W.W(out1,u,p,t)
  W.Z != nothing && W.Z(out2,u,p,t)
end

function NoiseFunction(t0,W,Z=nothing;kwargs...)
  iip=DiffEqBase.isinplace(W,4)
  NoiseFunction{iip}(t0,W,Z;kwargs...)
end

mutable struct NoiseGrid{T,N,Tt,T2,T3,ZType,inplace} <: AbstractNoiseProcess{T,N,Vector{T2},inplace}
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
(W::NoiseGrid)(u,p,t) = interpolate!(W,t)
(W::NoiseGrid)(out1,out2,u,p,t) = interpolate!(out1,out2,W,t)

mutable struct NoiseApproximation{T,N,Tt,T2,T3,S1,S2,ZType,inplace} <: AbstractNoiseProcess{T,N,Vector{T2},inplace}
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
(W::NoiseApproximation)(u,p,t) = interpolate!(W,t)
(W::NoiseApproximation)(out1,out2,u,p,t) = interpolate!(out1,out2,W,t)
