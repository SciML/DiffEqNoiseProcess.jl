DiffEqBase.has_reinit(i::AbstractNoiseProcess) = true
function DiffEqBase.reinit!(W::Union{NoiseProcess, NoiseApproximation}, dt;
                            t0 = W.t[1],
                            erase_sol = true,
                            setup_next = false)
    if erase_sol
        resize!(W.t, 1)
        resize!(W.W, 1)
        if W.Z != nothing
            resize!(W.Z, 1)
        end
    end

    W.curt = t0
    W.dt = dt
    if typeof(W) <: NoiseApproximation
        reinit!(W.source1)
        if W.source2 != nothing
            reinit!(W.source2)
        end
    end

    if isinplace(W)
        W.curW .= first(W.W)
        if W.Z != nothing
            W.curZ .= first(W.Z)
        end
    else
        W.curW = first(W.W)
        if W.Z != nothing
            W.curZ = first(W.Z)
        end
    end

    if typeof(W) <: NoiseProcess
        while !isempty(W.S₁)
            pop!(W.S₁) # Get a reset for this stack?
        end
        ResettableStacks.reset!(W.S₂)
    end
    setup_next && setup_next_step!(W)
    return nothing
end

function DiffEqBase.reinit!(W::AbstractNoiseProcess, dt;
                            t0 = W.t[1],
                            erase_sol = true,
                            setup_next = false)
    W.curt = t0
    W.dt = dt
    if typeof(W) <: NoiseGrid
        if t0 == W.t[1]
            idx = 1
        else
            tdir = sign(W.t[2] - W.t[1])
            @inbounds idx = searchsortedfirst(W.t, t0 - tdir * 10eps(typeof(t0)),
                                              rev = tdir < 0)
            # (tdir < 0) && (idx -= one(idx))
        end

        if isinplace(W)
            W.curW .= W.W[idx]
            if W.Z != nothing
                W.curZ .= W.Z[idx]
            end
        else
            W.curW = W.W[idx]
            if W.Z != nothing
                W.curZ = W.Z[idx]
            end
        end
        W.step_setup = true
    end
    setup_next && setup_next_step!(W)
    return nothing
end

function DiffEqBase.reinit!(W::NoiseFunction, dt;
                            t0 = W.t0,
                            erase_sol = true,
                            setup_next = false)
    W.curt = t0
    W.dt = dt

    if isinplace(W)
        W.curW .= first(W(t0))
        W.dW .= W.curW
        if W.Z !== nothing
            W.curZ .= last(W(t0))
            W.dZ .= W.curZ
        end
    else
        W.curW = first(W(t0))
        W.dW = W.curW
        if W.Z !== nothing
            W.curZ = last(W(t0))
            W.dZ = W.curZ
        end
    end

    return nothing
end

function Base.reverse(W::AbstractNoiseProcess)
    if typeof(W) <: NoiseGrid
        if W.Z === nothing
            backwardnoise = NoiseGrid(reverse(W.t), reverse(W.W))
        else
            backwardnoise = NoiseGrid(reverse(W.t), reverse(W.W), reverse(W.Z))
        end
    else
        W.save_everystep = false
        backwardnoise = NoiseWrapper(W, reverse = true)
    end
    return backwardnoise
end

function Base.similar(np::NoiseProcess, ::Type{NoiseProcess} = NoiseProcess)
    NoiseProcess{isinplace(np)}(0.0, 0.0, np.Z isa AbstractVector ? np.Z[1] : np.Z, np.dist,
                                np.bridge;
                                rswm = np.rswm, save_everystep = np.save_everystep,
                                rng = deepcopy(np.rng),
                                reset = np.reset, reseed = np.reseed,
                                continuous = np.continuous,
                                cache = np.cache)
end

function Base.copy(np::NoiseProcess)
    np2 = similar(np)
    for f in propertynames(np)
        setfield!(np2, f, getfield(np, f))
    end
    np2
end

function SciMLBase.remake(np::NoiseProcess; kwargs...)
    np_new = copy(np)
    inits = (t0 = :t, W0 = :W, Z0 = :Z)
    for kwarg in kwargs
        if first(kwarg) in keys(inits)
            setfield!(np_new, inits[first(kwarg)], [second(kwarg)])
        else
            setfield!(np_new, kwarg...)
        end
    end
    np_new
end
