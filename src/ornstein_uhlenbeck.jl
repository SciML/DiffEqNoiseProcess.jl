
struct OrnsteinUhlenbeck{T1,T2,T3}
  Θ::T1
  μ::T2
  σ::T3
end
# http://www.math.ku.dk/~susanne/StatDiff/Overheads1b.pdf
function (p::OrnsteinUhlenbeck)(W,dt,rng) #dist
  if typeof(W.dW) <: AbstractArray
    rand_val = wiener_randn(rng,W.dW)
  else
    rand_val = wiener_randn(rng,typeof(W.dW))
  end
  drift = p.μ .+ (W[end] .- p.μ) .* exp.(-p.Θ*dt)
  diffusion = p.σ .* sqrt.((1 .- exp.(-2p.Θ*dt))./(2p.Θ))
  drift .+ rand_val .* diffusion .- W[end]
end

#=
http://www.tandfonline.com/doi/pdf/10.1080/14697688.2014.941913?needAccess=true
http://www.tandfonline.com/doi/full/10.1080/14697688.2014.941913?src=recsys
=#

#=
(qX + r) = Θ(μ-x) = Θμ - Θx
q = -Θ
r = Θμ
https://arxiv.org/pdf/1011.0067.pdf page 18
=#
function ou_bridge(ou,W,W0,Wh,q,h) end
function ou_bridge!(rand_vec,ou,W,W0,Wh,q,h) end

function OrnsteinUhlenbeckProcess(Θ,μ,σ,t0,W0,Z0=nothing;kwargs...)
  ou = OrnsteinUhlenbeck(Θ,μ,σ)
  NoiseProcess(t0,W0,Z0,ou,nothing;kwargs...)
end

struct OrnsteinUhlenbeck!{T1,T2,T3}
  Θ::T1
  μ::T2
  σ::T3
end

function (p::OrnsteinUhlenbeck!)(rand_vec,W,dt,rng) #dist!
  wiener_randn!(rng,rand_vec)
  @. rand_vec = p.μ+(W[end]-p.μ)*exp(-p.Θ*dt) + rand_vec*p.σ*sqrt((1-exp.(-2*p.Θ.*dt))/(2*p.Θ)) - W[end]
end
function OrnsteinUhlenbeckProcess!(Θ,μ,σ,t0,W0,Z0=nothing;kwargs...)
  ou = OrnsteinUhlenbeck!(Θ,μ,σ)
  NoiseProcess(t0,W0,Z0,ou,nothing;kwargs...)
end
