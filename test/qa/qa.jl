using SciMLTesting, DiffEqNoiseProcess, JET, Test

run_qa(
    DiffEqNoiseProcess;
    explicit_imports = true,
    # ambiguities / deps_compat / piracies are genuine findings tracked in
    # https://github.com/SciML/DiffEqNoiseProcess.jl/issues/283
    aqua_broken = (:ambiguities, :deps_compat, :piracies),
    ei_kwargs = (;
        # AbstractNoiseProcess/AbstractNoiseProblem/DEIntegrator and @.. are
        # re-exported by DiffEqBase but owned by SciMLBase/FastBroadcast; they are
        # imported via DiffEqBase by long-standing convention (and are not declared
        # public in the owner module yet).
        all_explicit_imports_via_owners = (;
            ignore = (Symbol("@.."), :AbstractNoiseProblem, :AbstractNoiseProcess, :DEIntegrator),
        ),
        all_explicit_imports_are_public = (;
            ignore = (
                Symbol("@.."), :AbstractNoiseProblem, :AbstractNoiseProcess, :DEIntegrator,
                :step!,
            ),  # step! owned by CommonSolve, not declared public there
        ),
        # __solve / has_reinit are SciMLBase names re-exported by DiffEqBase and
        # extended here via DiffEqBase.<name>(...).
        all_qualified_accesses_via_owners = (; ignore = (:__solve, :has_reinit)),
        all_qualified_accesses_are_public = (;
            ignore = (
                Symbol("@pure"),        # Base
                :Broadcasted,           # Base.Broadcast
                :Experimental,          # Base
                :register_error_hint,   # Base.Experimental
                :AbstractDEAlgorithm,   # SciMLBase
                :ODE_DEFAULT_NORM,      # DiffEqBase
                :__solve,               # SciMLBase (extended via DiffEqBase)
                :has_reinit,            # SciMLBase (extended via DiffEqBase)
                :copyat_or_push!,       # ResettableStacks
                :default_rng, :seed!,   # Random — public on 1.11+, only flagged on the 1.10 LTS
                :sample_box, :sample_tail, :sample_wedge,  # own self-qualified names
            ),
        ),
    ),
)
