using DiffEqNoiseProcess
using Test

@time begin
  @time @testset "Interpolation Test" begin include("interpolation_test.jl") end
  @time @testset "RSwM1 Test" begin include("RSwM1_test.jl") end
  @time @testset "RSwM2 Test" begin include("RSwM2_test.jl") end
  @time @testset "RSwM3 Test" begin include("RSwM3_test.jl") end
  @time @testset "NoiseWrapper Test" begin include("noise_wrapper.jl") end
  @time @testset "NoiseFunction Test" begin include("noise_function.jl") end
  @time @testset "NoiseGrid Test" begin include("noise_grid.jl") end
  @time @testset "NoiseApproximation Test" begin include("noise_approximation.jl") end
  @time @testset "SDE NoiseWrapper Test" begin include("sde_noise_wrapper.jl") end
  @time @testset "Multidim Test" begin include("multi_dim.jl") end
  @time @testset "GBM Test" begin include("geometric_bm.jl") end
  @time @testset "OU Test" begin include("ornstein.jl") end
  @time @testset "Bridge Test" begin include("bridge_test.jl") end
  @time @testset "Adaptive SDE Distribution Test" begin include("sde_adaptivedistribution_tests.jl") end
end
