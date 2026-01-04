DiffEqBase.has_reinit(i::AbstractNoiseProcess) = true
function DiffEqBase.reinit!(
        W::Union{NoiseProcess, NoiseApproximation}, dt;
        t0 = W.t[1],
        erase_sol = true,
        setup_next = false
    )
    if erase_sol
        resize!(W.t, 1)
        resize!(W.W, 1)
        if W.Z !== nothing
            resize!(W.Z, 1)
        end
    end

    if W isa NoiseApproximation
        reinit!(W.source1)
        if W.source2 !== nothing
            reinit!(W.source2)
        end
    elseif W isa NoiseProcess
        ResettableStacks.reset!(W.S₁)
        while length(W.S₁) < length(W.reinitS₁)
            push!(W.S₁, W.reinitS₁.data[W.S₁.cur + 1])
        end
        ResettableStacks.reset!(W.S₂)
    end

    # Back to noise' starting state
    W.curt = first(W.t)

    if isinplace(W)
        W.curW .= first(W.W)
        if W.Z !== nothing
            W.curZ .= first(W.Z)
        end
    else
        W.curW = first(W.W)
        if W.Z !== nothing
            W.curZ = first(W.Z)
        end
    end

    if W.curt != t0
        # jump to prob's initial state if different from noise's
        W.dt = t0 - W.curt
        setup_next_step!(W, nothing, nothing)
        accept_step!(W, dt, nothing, nothing, false)
    end

    # prepare for actual first step
    W.dt = dt

    setup_next && setup_next_step!(W, nothing, nothing)

    return nothing
end

function DiffEqBase.reinit!(
        W::VirtualBrownianTree, dt;
        t0 = W.t[1],
        erase_sol = false,
        setup_next = true
    )

    # Back to noise's starting state
    W.curt = first(W.t)

    if isinplace(W)
        W.curW .= first(W.W)
        if W.Z !== nothing
            W.curZ .= first(W.Z)
        end
    else
        W.curW = first(W.W)
        if W.Z !== nothing
            W.curZ = first(W.Z)
        end
    end

    if W.curt != t0
        # jump to prob's initial state if different from noise's
        W.dt = t0 - W.curt
        setup_next_step!(W, nothing, nothing)
        accept_step!(W, dt, nothing, nothing, false)
    end

    # prepare for actual first step
    W.dt = dt
    setup_next && setup_next_step!(W, nothing, nothing)
    return nothing
end

function DiffEqBase.reinit!(
        W::NoiseGrid, dt;
        t0 = W.t[1],
        erase_sol = true,
        setup_next = false
    )
    W.curt = t0
    W.dt = dt
    if t0 == W.t[1]
        idx = 1
    else
        tdir = sign(W.t[2] - W.t[1])
        @inbounds idx = searchsortedfirst(
            W.t, t0 - tdir * 10eps(typeof(t0)),
            rev = tdir < 0
        )
        # (tdir < 0) && (idx -= one(idx))
    end

    if isinplace(W)
        W.curW .= W.W[idx]
        if W.Z !== nothing
            W.curZ .= W.Z[idx]
        end
    else
        W.curW = W.W[idx]
        if W.Z !== nothing
            W.curZ = W.Z[idx]
        end
    end
    W.step_setup = true
    setup_next && setup_next_step!(W, nothing, nothing)
    return nothing
end

function DiffEqBase.reinit!(
        W::AbstractNoiseProcess, dt;
        t0 = W.t[1],
        erase_sol = true,
        setup_next = false
    )
    W.curt = t0
    W.dt = dt

    setup_next && setup_next_step!(W, nothing, nothing)
    return nothing
end

function DiffEqBase.reinit!(
        W::NoiseFunction, dt;
        t0 = W.t0,
        erase_sol = true,
        setup_next = false
    )
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

function DiffEqBase.reinit!(
        W::NoiseTransport, dt;
        t0 = W.t0,
        erase_sol = true,
        setup_next = false
    )
    W.curt = t0
    W.dt = dt

    if W.rv isa AbstractArray
        W.RV(W.rng, W.rv)
    else
        W.rv = W.RV(W.rng)
    end

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

function DiffEqBase.reinit!(
        W::NoiseWrapper, dt;
        t0 = W.t[1],
        erase_sol = true,
        setup_next = false
    )
    W.curt = t0
    W.dt = dt

    if isinplace(W)
        W.curW .= W.W[1]
        W.dW .= W.curW
        if W.Z !== nothing
            W.curZ .= W.Z[1]
            W.dZ .= W.curZ
        end
    else
        W.curW = W.W[1]
        W.dW = W.curW
        if W.Z !== nothing
            W.curZ = W.Z[1]
            W.dZ = W.curZ
        end
    end

    return nothing
end

function Base.reverse(W::AbstractNoiseProcess)
    if W isa NoiseGrid
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
Base.reverse(W::VirtualBrownianTree) = W
