using DiffEqNoiseProcess


WHITE_NOISE_DIST  = (W,dt) -> sqrt(dt)*randn(size(W.dW))
WHITE_NOISE_BRIDGE= (W,W0,Wh,q,h) -> sqrt((1-q)*q*h)*randn(size(W.dW))+q*(Wh-W0)+W0
W = NoiseProcess(0.0,rand(4,4),nothing,WHITE_NOISE_DIST,WHITE_NOISE_BRIDGE,rswm=RSWM(adaptivealg=:RSwM3))

srand(200)
dt = 0.2
calculate_step!(W,dt)

for i in 1:100
  reject_step!(W,0.01)
  accept_step!(W,0.2)
  accept_step!(W,0.1)
end

dt = dt/100
for t in dt:dt:W.t[end]-dt
  W(t)
end

WHITE_NOISE_DIST  = function (rand_vec,W,dt)
  randn!(rand_vec)
  rand_vec .*= sqrt(dt)
end
WHITE_NOISE_BRIDGE = function (rand_vec,W,W0,Wh,q,h)
  randn!(rand_vec)
  rand_vec .= sqrt((1.-q).*q.*h).*rand_vec.+q.*(Wh.-W0).+W0
end
W = NoiseProcess(0.0,rand(4,4),rand(4,4),WHITE_NOISE_DIST,WHITE_NOISE_BRIDGE,rswm=RSWM(adaptivealg=:RSwM3))

srand(200)
dt = 0.2
calculate_step!(W,dt)

for i in 1:100
  reject_step!(W,0.01)
  accept_step!(W,0.2)
  accept_step!(W,0.1)
end

dt = dt/100
for t in dt:dt:W.t[end]-dt
  W(t)
end
