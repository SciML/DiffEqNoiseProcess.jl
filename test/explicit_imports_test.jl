using ExplicitImports
using DiffEqNoiseProcess
using Test

@testset "ExplicitImports" begin
    @test check_no_implicit_imports(DiffEqNoiseProcess) === nothing
    @test check_no_stale_explicit_imports(DiffEqNoiseProcess) === nothing
end
