using SciMLTesting, DiffEqNoiseProcess, JET, Test

function _public_marker_name(token)
    token = strip(token)
    token = replace(token, r"^\s*:" => "")
    token = replace(token, r"^Symbol\(\"(.+)\"\)$" => s"\1")
    return Symbol(token)
end

function _parse_public_marker_names(names)
    return [_public_marker_name(token) for token in split(names, ",") if !isempty(strip(token))]
end

function _declared_public_api()
    api = Set{Symbol}()
    src_root = joinpath(pkgdir(DiffEqNoiseProcess), "src")
    for (root, _, files) in walkdir(src_root)
        for file in files
            endswith(file, ".jl") || continue
            for line in eachline(joinpath(root, file))
                line = strip(split(line, '#'; limit = 2)[1])
                for marker in ("export", "public", "@public")
                    m = match(Regex("^" * marker * "\\s+(.+)\$"), line)
                    m === nothing || union!(api, _parse_public_marker_names(m.captures[1]))
                end
                for m in eachmatch(r"Expr\(\s*:public\s*,\s*(.*?)\)", line)
                    union!(api, _parse_public_marker_names(m.captures[1]))
                end
            end
        end
    end
    return sort!(collect(api))
end

function _rendered_docs_entries()
    entries = Set{Symbol}()
    docs_root = joinpath(pkgdir(DiffEqNoiseProcess), "docs", "src")
    for (root, _, files) in walkdir(docs_root)
        for file in files
            endswith(file, ".md") || continue
            in_docs_block = false
            for line in eachline(joinpath(root, file))
                line = strip(line)
                if line == "```@docs"
                    in_docs_block = true
                    continue
                elseif startswith(line, "```")
                    in_docs_block = false
                    continue
                end
                in_docs_block || continue
                isempty(line) && continue
                name = replace(line, r"\(.*" => "")
                name = split(name, ".")[end]
                push!(entries, Symbol(name))
            end
        end
    end
    return entries
end

@testset "public API documentation coverage" begin
    public_api = _declared_public_api()

    missing_docstrings = filter(name -> !Docs.hasdoc(DiffEqNoiseProcess, name), public_api)
    @test missing_docstrings == Symbol[]

    rendered_entries = _rendered_docs_entries()
    missing_rendered_entries = filter(name -> name ∉ rendered_entries, public_api)
    @test missing_rendered_entries == Symbol[]
end

run_qa(
    DiffEqNoiseProcess;
    explicit_imports = true,
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
