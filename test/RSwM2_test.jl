@testset "RSwM2" begin

using DiffEqNoiseProcess

W = WienerProcess(0.0,0.0,rswm=RSWM(adaptivealg=:RSwM2))

dt = 0.2
calculate_step!(W,dt,nothing,nothing)
W.dt,W.curt,W.curW,W.dW

reject_step!(W,0.1,nothing,nothing)
W.dt,W.curt,W.curW,W.dW
DiffEqNoiseProcess.first(W.S₁)

accept_step!(W,0.2,nothing,nothing)
W.dt,W.curt,W.curW,W.dW

reject_step!(W,0.1,nothing,nothing)
W.dt,W.curt,W.curW,W.dW
DiffEqNoiseProcess.first(W.S₁)

for i in 1:100
  reject_step!(W,0.05,nothing,nothing)
  accept_step!(W,0.1,nothing,nothing)
end

dt = dt/100
for t in dt:dt:W.t[end]-dt
  W(t)
end

end
