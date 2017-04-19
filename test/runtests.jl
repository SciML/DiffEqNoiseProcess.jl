using DiffEqNoiseProcess
using Base.Test

@time @testset "Interpolation Test" begin include("interpolation_test.jl") end
@time @testset "RSwM1 Test" begin include("RSwM1_test.jl") end
@time @testset "RSwM3 Test" begin include("RSwM3_test.jl") end
