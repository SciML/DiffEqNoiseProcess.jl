@testset "NoiseWrapper" begin

using DiffEqNoiseProcess, Test, Random

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


_W = WienerProcess(0.0,0.0,0.0)

dt = 0.1
calculate_step!(_W,dt)

for i in 1:10
  accept_step!(_W,dt)
end

old_W = copy(_W.W)

W2 = NoiseWrapper(_W)
lspace =range(_W.t[1], stop=_W.t[end], length=1000)
dt = lspace[2]-lspace[1]
calculate_step!(W2,dt)
for t in lspace
  accept_step!(W2,dt)
end

@test W2.W[end] ≈ _W(W2.t[end])[1]
@test W2.Z[end] ≈ _W(W2.t[end])[2]


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


_W = WienerProcess!(0.0,zeros(4),zeros(4))

dt = 0.1
calculate_step!(_W,dt)

for i in 1:10
  accept_step!(_W,dt)
end

W2 = NoiseWrapper(_W)
lspace =range(_W.t[1], stop=_W.t[end], length=20)
dt = lspace[2]-lspace[1]
calculate_step!(W2,dt)
for t in lspace
  accept_step!(W2,dt)
end

@test W2.W[end]≈ _W(W2.t[end])[1]
@test W2.Z[end]≈ _W(W2.t[end])[2]

@test W2.W[end] ≈ _W.W[end-1]

end
