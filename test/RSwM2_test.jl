@testset "RSwM2" begin

using DiffEqNoiseProcess

W = WienerProcess(0.0,0.0,rswm=RSWM(adaptivealg=:RSwM2))

dt = 0.2
calculate_step!(W,dt)
W.dt,W.curt,W.curW,W.dW

reject_step!(W,0.1)
W.dt,W.curt,W.curW,W.dW
DiffEqNoiseProcess.top(W.S₁)

accept_step!(W,0.2)
W.dt,W.curt,W.curW,W.dW

reject_step!(W,0.1)
W.dt,W.curt,W.curW,W.dW
DiffEqNoiseProcess.top(W.S₁)

for i in 1:100
  reject_step!(W,0.05)
  accept_step!(W,0.1)
end

dt = dt/100
for t in dt:dt:W.t[end]-dt
  W(t)
end

end
