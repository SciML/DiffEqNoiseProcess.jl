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
                Symbol("@.."),          # FastBroadcast (via DiffEqBase), not public
                :AbstractNoiseProblem,  # SciMLBase (via DiffEqBase), not public
                :AbstractNoiseProcess,  # SciMLBase (via DiffEqBase), not public
                :DEIntegrator,          # SciMLBase (via DiffEqBase), not public
            ),
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
                :ODE_DEFAULT_NORM,      # DiffEqBase, not public
                :__solve,               # DiffEqBase (SciMLBase name extended via DiffEqBase)
                :has_reinit,            # DiffEqBase (SciMLBase name extended via DiffEqBase)
                :copyat_or_push!,       # ResettableStacks, not public
            ),
        ),
    ),
)
