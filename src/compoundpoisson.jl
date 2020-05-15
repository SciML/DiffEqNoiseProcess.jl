struct CompoundPoissonProcess{R}
  rate::R
end
function (P::CompoundPoissonProcess)(W,dt,u,p,t,rng)
  PoissonRandom.pois_rand.(rng,dt*P.rate(u,p,t))
end

struct CompoundPoissonProcess!{R}
  rate::R
end
function (P::CompoundPoissonProcess!)(rand_vec,W,dt,u,p,t,rng)
  P.rate(rand_vec,u,p,t)
  @. rand_vec = PoissonRandom.pois_rand(rng,dt*rand_vec)
end

# https://www.math.wisc.edu/~anderson/papers/AndPostleap.pdf
# Incorporating postleap checks in tau-leaping
# J. Chem. Phys. 128, 054103 (2008); https://doi.org/10.1063/1.2819665
function cpp_bridge(cpp,W,W0,Wh,q,h,u,p,t,rng)
  rand.(rng,Distributions.Binomial.(Int.(Wh),float.(q)))
end
function cpp_bridge!(rand_vec,cpp,W,W0,Wh,q,h,u,p,t,rng)
  rand_vec .= rand.(rng,Distributions.Binomial.(Int.(Wh),float.(q)))
end

function CompoundPoissonProcess(rate,t0,W0,Z0=nothing;kwargs...)
  cpp = CompoundPoissonProcess(rate)
  NoiseProcess(t0,W0,Z0,cpp,(W,W0,Wh,q,h,u,p,t,rng)->cpp_bridge(cpp,W,W0,Wh,q,h,u,p,t,rng);continuous=false,kwargs...)
end

function CompoundPoissonProcess!(rate,t0,W0,Z0=nothing;kwargs...)
  cpp = CompoundPoissonProcess!(rate)
  NoiseProcess(t0,W0,Z0,cpp,(rand_vec,W,W0,Wh,q,h,u,p,t,rng)->cpp_bridge!(rand_vec,cpp,W,W0,Wh,q,h,u,p,t,rng);continuous=false,kwargs...)
end
