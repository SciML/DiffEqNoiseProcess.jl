struct GeometricBrownianMotion{T1,T2}
  μ::T1
  σ::T2
end

function (X::GeometricBrownianMotion)(dW,W,dt,u,p,t,rng) #dist
  drift = X.μ-(1/2)*X.σ^2
  if typeof(dW) <: AbstractArray
    rand_val = wiener_randn(rng,dW)
  else
    rand_val = wiener_randn(rng,typeof(dW))
  end
  new_val = @.. exp(drift*dt + X.σ*sqrt(dt)*rand_val)
  return W[end]*(new_val-1)
end

#=
u = q*h
t = h
s = 0

t-u = h-qh = (1-q)h
t-s = h
u - s = q*h
https://math.stackexchange.com/questions/412470/conditional-distribution-in-brownian-motion

=#
function gbm_bridge(dW,gbm,W,W0,Wh,q,h,u,p,t,rng)
  if typeof(dW) <: AbstractArray
    return gbm.σ*sqrt((1-q)*q*abs(h))*wiener_randn(rng,dW)+q*Wh
  else
    return gbm.σ*sqrt((1-q)*q*abs(h))*wiener_randn(rng,typeof(dW))+q*Wh
  end
end
function gbm_bridge!(rand_vec,gbm,W,W0,Wh,q,h,u,p,t,rng)
  wiener_randn!(rng,rand_vec)
  @.. rand_vec = gbm.σ*sqrt((1-q)*q*abs(h))*rand_vec+q*Wh
end

function GeometricBrownianMotionProcess(μ,σ,t0,W0,Z0=nothing;kwargs...)
  gbm = GeometricBrownianMotion(μ,σ)
  NoiseProcess{false}(t0,W0,Z0,gbm,(dW,W,W0,Wh,q,h,u,p,t,rng)->gbm_bridge(dW,gbm,W,W0,Wh,q,h,u,p,t,rng);kwargs...)
end

struct GeometricBrownianMotion!{T1,T2}
  μ::T1
  σ::T2
end
function (X::GeometricBrownianMotion!)(rand_vec,W,dt,u,p,t,rng) #dist!
  wiener_randn!(rng,rand_vec)
  @.. rand_vec = W[end]*(exp(X.μ-(1/2)*X.σ*dt + X.σ*sqrt(dt)*rand_vec)-1)
end
function GeometricBrownianMotionProcess!(μ,σ,t0,W0,Z0=nothing;kwargs...)
  gbm = GeometricBrownianMotion!(μ,σ)
  NoiseProcess{true}(t0,W0,Z0,gbm,(rand_vec,W,W0,Wh,q,h,u,p,t,rng)->gbm_bridge!(rand_vec,gbm,W,W0,Wh,q,h,u,p,t,rng);kwargs...)
end
