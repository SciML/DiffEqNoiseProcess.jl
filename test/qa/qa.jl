using SciMLTesting, DiffEqNoiseProcess, JET, Test

run_qa(
    DiffEqNoiseProcess;
    # ambiguities / deps_compat / piracies are genuine findings tracked in
    # https://github.com/SciML/DiffEqNoiseProcess.jl/issues/283
    aqua_broken = (:ambiguities, :deps_compat, :piracies),
    # On Julia >= 1.12, JET typo-mode (report_package) flags 27 pre-existing
    # "local variable may be undefined" reports, all of the same shape: a value
    # assigned inside one `X !== nothing` / `offset !== nothing` guard and used
    # inside a *second* branch protected by the identical (never-invalidated)
    # guard, in interpolate! / generate_boxes / the VBT bridge. JET's
    # intraprocedural undef analysis can't prove the guards are correlated, so
    # they are false positives at runtime but genuine latent code patterns. The
    # previous hand-rolled qa explicitly skipped report_package (only report_opt
    # on a few functions), so these were never surfaced; they reproduce
    # byte-identically on unmodified master. On Julia <= 1.11 JET finds 0 reports
    # (the analysis is version-specific), so jet_broken is gated to avoid an
    # Unexpected Pass there.
    # Tracked in https://github.com/SciML/DiffEqNoiseProcess.jl/issues/283
    jet_broken = VERSION >= v"1.12",
    ei_kwargs = (;
        # @.. (FastBroadcast) and DEIntegrator (SciMLBase) are re-exported by
        # DiffEqBase but owned elsewhere; imported via DiffEqBase by convention.
        all_explicit_imports_via_owners = (;
            ignore = (Symbol("@.."), :DEIntegrator),
        ),
        all_explicit_imports_are_public = (;
            ignore = (
                Symbol("@.."),   # FastBroadcast (via DiffEqBase), not public
                :DEIntegrator,   # SciMLBase (via DiffEqBase), not public
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
                :__solve,               # DiffEqBase (SciMLBase name extended via DiffEqBase)
                :has_reinit,            # DiffEqBase (SciMLBase name extended via DiffEqBase)
                :copyat_or_push!,       # ResettableStacks, not public
            ),
        ),
    ),
)
