struct OrnsteinUhlenbeck{T1,T2,T3}
  Θ::T1
  μ::T2
  σ::T3
end
# http://www.math.ku.dk/~susanne/StatDiff/Overheads1b.pdf
function (X::OrnsteinUhlenbeck)(dW,W,dt,u,p,t,rng) #dist
  if typeof(dW) <: AbstractArray
    rand_val = wiener_randn(rng,dW)
  else
    rand_val = wiener_randn(rng,typeof(dW))
  end
  drift = X.μ .+ (W[end] .- X.μ) .* exp.(-X.Θ*dt)
  diffusion = X.σ .* sqrt.((1 .- exp.(-2X.Θ*dt))./(2X.Θ))
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
function ou_bridge(dW,ou,W,W0,Wh,q,hu,p,t,rng) end
function ou_bridge!(rand_vec,ou,W,W0,Wh,q,h,u,p,t,rng) end

function OrnsteinUhlenbeckProcess(Θ,μ,σ,t0,W0,Z0=nothing;kwargs...)
  ou = OrnsteinUhlenbeck(Θ,μ,σ)
  NoiseProcess{false}(t0,W0,Z0,ou,nothing;kwargs...)
end

struct OrnsteinUhlenbeck!{T1,T2,T3}
  Θ::T1
  μ::T2
  σ::T3
end

function (X::OrnsteinUhlenbeck!)(rand_vec,W,dt,u,p,t,rng) #dist!
  wiener_randn!(rng,rand_vec)
  @.. rand_vec = X.μ+(W[end]-X.μ)*exp(-X.Θ*dt) + rand_vec*X.σ*sqrt((1-exp.(-2*X.Θ.*dt))/(2*X.Θ)) - W[end]
end
function OrnsteinUhlenbeckProcess!(Θ,μ,σ,t0,W0,Z0=nothing;kwargs...)
  ou = OrnsteinUhlenbeck!(Θ,μ,σ)
  NoiseProcess{true}(t0,W0,Z0,ou,nothing;kwargs...)
end
