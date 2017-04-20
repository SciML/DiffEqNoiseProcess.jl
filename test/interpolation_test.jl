using DiffEqNoiseProcess

WHITE_NOISE_DIST  = (W,dt) -> sqrt(dt)*randn()
WHITE_NOISE_BRIDGE= (W,W0,Wh,q,h) -> sqrt((1-q)*q*h)*randn()+q*(Wh-W0)+W0
W = NoiseProcess(0.0,0.0,0.0,WHITE_NOISE_DIST,WHITE_NOISE_BRIDGE,rswm=RSWM(adaptivealg=:RSwM1))

dt = 0.1
calculate_step!(W,dt)

for i in 1:10
  accept_step!(W,dt)
end

dt = dt/100
for t in dt:dt:1-dt
  W(t)
end

dt = dt/100
for t in dt:dt:1-dt
  W(t)
end
