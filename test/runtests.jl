using Test, Pkg

const GROUP = get(ENV, "GROUP", "All")

function activate_gpu_env()
    Pkg.activate("gpu")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    return Pkg.instantiate()
end

function activate_qa_env()
    Pkg.activate("qa")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    return Pkg.instantiate()
end

@time begin
    if GROUP == "All" || GROUP == "Core1"
        include("explicit_imports_test.jl")
        include("interpolation_test.jl")
        include("RSwM1_test.jl")
        include("RSwM2_test.jl")
        include("RSwM3_test.jl")
        include("correlated.jl")
        include("noise_wrapper.jl")
        include("noise_function.jl")
        include("noise_transport.jl")
        include("copy_noise_test.jl")
        include("VBT_test.jl")
        include("noise_grid.jl")
        include("noise_approximation.jl")
        include("sde_noise_wrapper.jl")
        include("multi_dim.jl")
        include("geometric_bm.jl")
        include("compoundpoisson.jl")
        include("ornstein.jl")
        include("ensemble_test.jl")
        include("sde_adaptivedistribution_tests.jl")
        include("reversal_test.jl")
        include("two_processes.jl")
        include("extraction_test.jl")
        include("restart_test.jl")
        include("reinit_test.jl")
        include("BWT_test.jl")
        include("pcn_test.jl")
        include("savestep_test.jl")
        include("simple_noise_bigfloat_test.jl")
    end

    if GROUP == "All" || GROUP == "Bridge"
        include("bridge_test.jl")
    end

    if GROUP == "QA"
        activate_qa_env()
        include("qa/qa_tests.jl")
    end
end
