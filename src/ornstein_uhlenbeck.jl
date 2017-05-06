function ou_dist(W,dt)
  if typeof(W.dW) <: AbstractArray
    return sqrt(abs(dt))*wiener_randn(size(W.dW))
  else
    return sqrt(abs(dt))*wiener_randn(typeof(W.dW))
  end
end
#=
function ou_bridge(W,W0,Wh,q,h)
  sqrt((1-q)*q*abs(h))*wiener_randn(typeof(W.dW))+q*(Wh-W0)+W0
end
=#
OrnsteinUhlenbeckProcess(t0,W0,Z0=nothing) = NoiseProcess(t0,W0,Z0,ou_dist,nothing,rswm=RSWM())

function ou_dist!(rand_vec,W,dt)
  wiener_randn!(rand_vec)
  rand_vec .*= sqrt(abs(dt))
end
#=
function ou_bridge!(rand_vec,W,W0,Wh,q,h)
  wiener_randn!(rand_vec)
  rand_vec .= sqrt((1.-q).*q.*abs(h)).*rand_vec.+q.*(Wh.-W0).+W0
end
=#
OrnsteinUhlenbeckProcess!(t0,W0,Z0=nothing) = NoiseProcess(t0,W0,Z0,ou_dist!,nothing,rswm=RSWM())
