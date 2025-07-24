using DiffEqNoiseProcess
using Documenter

# Test if we can find the docstrings
println("Testing docstring availability...")

# Check exported functions
for fname in [:accept_step!, :reject_step!, :calculate_step!, :setup_next_step!, :save_noise!]
    if isdefined(DiffEqNoiseProcess, fname)
        println("✓ $fname is defined")
        docs = Docs.doc(getfield(DiffEqNoiseProcess, fname))
        if !isempty(string(docs))
            println("  ✓ Has documentation")
        else
            println("  ✗ No documentation found")
        end
    else
        println("✗ $fname is not defined")
    end
end

# Check types
println("\nChecking type documentation...")
for tname in [:RSWM, :WienerProcess, :NoiseProcess]
    if isdefined(DiffEqNoiseProcess, tname)
        println("✓ $tname is defined")
        docs = Docs.doc(getfield(DiffEqNoiseProcess, tname))
        if !isempty(string(docs)) && !occursin("No documentation found", string(docs))
            println("  ✓ Has documentation")
        else
            println("  ✗ No documentation found")
        end
    else
        println("✗ $tname is not defined")
    end
end