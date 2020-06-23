# https://www.math.wisc.edu/~anderson/papers/AndPostleap.pdf
# Incorporating postleap checks in tau-leaping
# J. Chem. Phys. 128, 054103 (2008); https://doi.org/10.1063/1.2819665
function cpp_bridge(dW,cpp,W,W0,Wh,q,h,u,p,t,rng)
  rand.(rng,Distributions.Binomial.(Int.(Wh),float.(q)))
end
function cpp_bridge!(rand_vec,cpp,W,W0,Wh,q,h,u,p,t,rng)
  rand_vec .= rand.(rng,Distributions.Binomial.(Int.(Wh),float.(q)))
end

mutable struct CompoundPoissonProcess{R,CR}
  rate::R
  currate::CR
  computerates::Bool
  function CompoundPoissonProcess(rate,t0,W0;computerates=true,kwargs...)
    cpp = new{typeof(rate),typeof(W0)}(rate,W0,computerates)
    NoiseProcess{false}(t0,W0,nothing,cpp,(dW,W,W0,Wh,q,h,u,p,t,rng)->cpp_bridge(dW,cpp,W,W0,Wh,q,h,u,p,t,rng);continuous=false,cache=cpp,kwargs...)
  end
end
function (P::CompoundPoissonProcess)(dW,W,dt,u,p,t,rng)
  P.computerates && (P.currate = P.rate(u,p,t))
  PoissonRandom.pois_rand.(rng,dt.*P.currate)
end

struct CompoundPoissonProcess!{R,CR}
  rate::R
  currate::CR
  computerates::Bool
  function CompoundPoissonProcess!(rate,t0,W0;computerates=true,kwargs...)
    cpp = new{typeof(rate),typeof(W0)}(rate,copy(W0),computerates)
    NoiseProcess{true}(t0,W0,nothing,cpp,(rand_vec,W,W0,Wh,q,h,u,p,t,rng)->cpp_bridge!(rand_vec,cpp,W,W0,Wh,q,h,u,p,t,rng);continuous=false,cache=cpp,kwargs...)
  end
end
function (P::CompoundPoissonProcess!)(rand_vec,W,dt,u,p,t,rng)
  P.computerates && P.rate(P.currate,u,p,t)
  @. rand_vec = PoissonRandom.pois_rand(rng,dt*P.currate)
end
