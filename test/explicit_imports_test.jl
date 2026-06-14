@safetestset "ExplicitImports" begin
    using ExplicitImports
    using DiffEqNoiseProcess
    using Test
    @test check_no_implicit_imports(DiffEqNoiseProcess) === nothing
    @test check_no_stale_explicit_imports(DiffEqNoiseProcess) === nothing
end
