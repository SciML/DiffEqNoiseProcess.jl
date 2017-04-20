using DiffEqNoiseProcess
using Base.Test

@time @testset "Interpolation Test" begin include("interpolation_test.jl") end
@time @testset "RSwM1 Test" begin include("RSwM1_test.jl") end
@time @testset "RSwM2 Test" begin include("RSwM2_test.jl") end
@time @testset "RSwM3 Test" begin include("RSwM3_test.jl") end
@time @testset "NoiseWrapper Test" begin include("noise_wrapper.jl") end
@time @testset "SDE NoiseWrapper Test" begin include("sde_noise_wrapper.jl") end
@time @testset "Multidim Test" begin include("multi_dim.jl") end
@time @testset "Adaptive SDE Distribution Test" begin include("sde_adaptivedistribution_tests.jl") end
