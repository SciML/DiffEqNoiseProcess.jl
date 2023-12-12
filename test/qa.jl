using DiffEqNoiseProcess, Aqua
@testset "Aqua" begin
    Aqua.find_persistent_tasks_deps(DiffEqNoiseProcess)
    Aqua.test_ambiguities(DiffEqNoiseProcess, recursive = false)
    Aqua.test_deps_compat(DiffEqNoiseProcess)
    Aqua.test_piracies(DiffEqNoiseProcess,
        treat_as_own = [DiffEqNoiseProcess.AbstractNoiseProcess,
            DiffEqNoiseProcess.AbstractDEAlgorithm])
    Aqua.test_project_extras(DiffEqNoiseProcess)
    Aqua.test_stale_deps(DiffEqNoiseProcess)
    Aqua.test_unbound_args(DiffEqNoiseProcess)
    Aqua.test_undefined_exports(DiffEqNoiseProcess)
end
