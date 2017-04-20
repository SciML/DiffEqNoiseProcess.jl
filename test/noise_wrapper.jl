using DiffEqNoiseProcess, Base.Test

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

W2 = NoiseWrapper(_W)
lspace =linspace(_W.t[1],_W.t[end],10000)
dt = lspace[2]-lspace[1]
calculate_step!(W2,dt)
for t in lspace
  accept_step!(W2,dt)
end

@test W2.W[end]≈ _W(W2.t[end])[1]
@test W2.Z[end]≈ _W(W2.t[end])[2]


# Inplace


_W = WienerProcess!(0.0,zeros(4),zeros(4))

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

W2 = NoiseWrapper(_W)
lspace =linspace(_W.t[1],_W.t[end],10000)
dt = lspace[2]-lspace[1]
calculate_step!(W2,dt)
for t in lspace
  accept_step!(W2,dt)
end

@test W2.W[end]≈ _W(W2.t[end])[1]
@test W2.Z[end]≈ _W(W2.t[end])[2]

@test W2.W[end] ≈ _W.W[end-1]
