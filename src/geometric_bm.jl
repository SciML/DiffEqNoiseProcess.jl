immutable GeometricBrownianMotion{T1,T2}
  μ::T1
  σ::T2
end

function (p::GeometricBrownianMotion)(W,dt) #dist
  drift = p.μ-(1/2)*p.σ^2
  if typeof(W.dW) <: AbstractArray
    rand_val = wiener_randn(size(W.dW))
  else
    rand_val = wiener_randn(typeof(W.dW))
  end
  new_val = exp(drift*dt + p.σ*sqrt(dt)*rand_val)
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
function gbm_bridge(gbm,W,W0,Wh,q,h)
  if typeof(W.dW) <: AbstractArray
    return gbm.σ*sqrt((1-q)*q*abs(h))*wiener_randn(size(W.dW))+q*(Wh-W0)+W0
  else
    return gbm.σ*sqrt((1-q)*q*abs(h))*wiener_randn(typeof(W.dW))+q*(Wh-W0)+W0
  end
end
function gbm_bridge!(rand_vec,gbm,W,W0,Wh,q,h)
  wiener_randn!(rand_vec)
  rand_vec .= gbm.σ.*sqrt((1.-q).*q.*abs(h)).*rand_vec.+q.*(Wh.-W0).+W0
end

function GeometricBrownianMotionProcess(μ,σ,t0,W0,Z0=nothing)
  gbm = GeometricBrownianMotion(μ,σ)
  NoiseProcess(t0,W0,Z0,gbm,(W,W0,Wh,q,h)->gbm_bridge(gbm,W,W0,Wh,q,h),rswm=RSWM())
end

immutable GeometricBrownianMotion!{T1,T2}
  μ::T1
  σ::T2
end
function (p::GeometricBrownianMotion!)(rand_vec,W,dt) #dist!
  wiener_randn!(rand_vec)
  rand_vec .= W[end].*(exp.(p.μ.-(1/2).*p.σ.*dt .+ p.σ.*sqrt(dt).*rand_vec).-1)
end
function GeometricBrownianMotionProcess!(μ,σ,t0,W0,Z0=nothing)
  gbm = GeometricBrownianMotion!(μ,σ)
  NoiseProcess(t0,W0,Z0,gbm,(rand_vec,W,W0,Wh,q,h)->gbm_bridge!(rand_vec,gbm,W,W0,Wh,q,h),rswm=RSWM())
end
