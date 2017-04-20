using DiffEqNoiseProcess

WHITE_NOISE_DIST  = (W,dt) -> sqrt(dt)*randn()
WHITE_NOISE_BRIDGE= (W,W0,Wh,q,h) -> sqrt((1-q)*q*h)*randn()+q*(Wh-W0)+W0
W = NoiseProcess(0.0,0.0,WHITE_NOISE_DIST,WHITE_NOISE_BRIDGE,rswm=RSWM(adaptivealg=:RSwM1))

dt = 0.2
calculate_step!(W,dt)
W.dt,W.curt,W.curW,W.dW

reject_step!(W,0.1)
W.dt,W.curt,W.curW,W.dW
DiffEqNoiseProcess.top(W.S₁)

accept_step!(W,0.2)
W.dt,W.curt,W.curW,W.dW

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
