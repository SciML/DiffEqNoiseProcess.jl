using DiffEqNoiseProcess

_W = WienerProcess(0.0,0.0,0.0)

dt = 0.1
calculate_step!(_W,dt)

for i in 1:10
  accept_step!(_W,dt)
end

W2 = NoiseWrapper(_W)

dt = 0.1
calculate_step!(W2,dt)

for i in 1:10
  accept_step!(W2,dt)
end

using Plots
plot(_W.t,_W.W)

W2 = NoiseWrapper(_W)
lspace =linspace(_W.t[1],_W.t[end],10000)
dt = lspace[2]-lspace[1]
calculate_step!(W2,dt)
for t in lspace
  accept_step!(W2,dt)
end

plot!(W2.t,W2.W)
