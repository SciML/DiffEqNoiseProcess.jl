@testset "Multidim" begin

using DiffEqNoiseProcess, Random, Statistics

W = WienerProcess(0.0,rand(4,4),rswm=RSWM(adaptivealg=:RSwM3))

Random.seed!(200)
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

W = WienerProcess!(0.0,rand(4,4),rand(4,4),rswm=RSWM(adaptivealg=:RSwM3))

Random.seed!(200)
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

end
