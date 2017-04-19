using DiffEqNoiseProcess

WHITE_NOISE_DIST  = (W,dt) -> sqrt(dt)*randn()
WHITE_NOISE_BRIDGE= (W0,Wh,q,h) -> sqrt((1-q)*q*h)*randn()+q*(Wh-W0)+W0
W = WienerProcess(0.0,0.0,WHITE_NOISE_DIST,WHITE_NOISE_BRIDGE,rswm=StochasticDiffEq.RSWM(adaptivealg=:RSwM3))

dt = 0.2
calculate_step!(W,dt)
W.curW,W.dW,W.curZ,W.dZ
isempty(W.S₂)

reject_step!(W,0.1)
W.curW,W.dW,W.curZ,W.dZ
DiffEqNoiseProcess.top(W.S₁)
isempty(W.S₂)

accept_step!(W,0.2)
W.curW,W.dW,W.curZ,W.dZ
isempty(W.S₂)

reject_step!(W,0.1)
W.curW,W.dW,W.curZ,W.dZ
DiffEqNoiseProcess.top(W.S₁)

accept_step!(W,0.2)
W.curW,W.dW,W.curZ,W.dZ

for i in 1:100
  reject_step!(W,0.1)
  accept_step!(W,0.1)
end

using Plots
plot(W.t,W.W)
plot!(W.t,W.Z)

dt = dt/100
for t in dt:dt:W.t[end]-dt
  W(t)
end

plot(W.t,W.W)
plot!(W.t,W.Z)





WHITE_NOISE_DIST  = (W,dt) -> sqrt(dt)*randn()
WHITE_NOISE_BRIDGE= (W0,Wh,q,h) -> sqrt((1-q)*q*h)*randn()+q*(Wh-W0)+W0
W = WienerProcess(0.0,0.0,WHITE_NOISE_DIST,WHITE_NOISE_BRIDGE,rswm=StochasticDiffEq.RSWM(adaptivealg=:RSwM3))

dt = 0.2
calculate_step!(W,dt)
W.curW,W.dW,W.curZ,W.dZ
isempty(W.S₂)

for i in 1:100
  reject_step!(W,0.01)
  accept_step!(W,0.2)
  accept_step!(W,0.1)
end

using Plots
plot(W.t,W.W)
plot!(W.t,W.Z)

dt = dt/100
for t in dt:dt:W.t[end]-dt
  W(t)
end

plot(W.t,W.W)
plot!(W.t,W.Z)
