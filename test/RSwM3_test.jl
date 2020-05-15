@testset "RSwM3" begin

using DiffEqNoiseProcess

W = WienerProcess(0.0,0.0,0.0,rswm=RSWM(adaptivealg=:RSwM3))

dt = 0.2
calculate_step!(W,dt,nothing,nothing)
W.curW,W.dW,W.curZ,W.dZ
isempty(W.S₂)

reject_step!(W,0.1,nothing,nothing)
W.curW,W.dW,W.curZ,W.dZ
DiffEqNoiseProcess.first(W.S₁)
isempty(W.S₂)

accept_step!(W,0.2,nothing,nothing)
W.curW,W.dW,W.curZ,W.dZ
isempty(W.S₂)

reject_step!(W,0.1,nothing,nothing)
W.curW,W.dW,W.curZ,W.dZ
DiffEqNoiseProcess.first(W.S₁)

accept_step!(W,0.2,nothing,nothing)
W.curW,W.dW,W.curZ,W.dZ

for i in 1:100
  reject_step!(W,0.1,nothing,nothing)
  accept_step!(W,0.1,nothing,nothing)
end

dt = dt/100
for t in dt:dt:W.t[end]-dt
  W(t)
end


W = WienerProcess(0.0,0.0,0.0,rswm=RSWM(adaptivealg=:RSwM3))

dt = 0.2
calculate_step!(W,dt,nothing,nothing)
W.curW,W.dW,W.curZ,W.dZ
isempty(W.S₂)

for i in 1:100
  reject_step!(W,0.01,nothing,nothing)
  accept_step!(W,0.2,nothing,nothing)
  accept_step!(W,0.1,nothing,nothing)
end

dt = dt/100
for t in dt:dt:W.t[end]-dt
  W(t)
end

end
