function Base.copy!(Wnew::T, W::T) where {T <: AbstractNoiseProcess}
    for x in filter(!=(:u), fieldnames(typeof(W)))
        if !ismutable(getfield(W, x))
            setfield!(Wnew, x, getfield(W, x))
        elseif getfield(W, x) isa AbstractNoiseProcess
            copy!(getfield(Wnew, x), getfield(W, x))
        elseif getfield(W, x) isa AbstractArray && !ismutable(eltype(getfield(W, x)))
            setfield!(Wnew, x, copy(getfield(W, x)))
        elseif getfield(W, x) isa AbstractArray
            setfield!(Wnew, x, recursivecopy(getfield(W, x)))
        elseif getfield(W, x) isa ResettableStacks.ResettableStack
            setfield!(getfield(Wnew, x), :cur, getfield(W, x).cur)
            setfield!(getfield(Wnew, x), :numResets, getfield(W, x).numResets)
            setfield!(getfield(Wnew, x), :data, recursivecopy(getfield(W, x).data))
        elseif getfield(W, x) isa RSWM
            setfield!(getfield(Wnew, x), :discard_length, getfield(W, x).discard_length)
            setfield!(getfield(Wnew, x), :adaptivealg, getfield(W, x).adaptivealg)
        elseif typeof(getfield(W, x)) <:
               Union{BoxGeneration1, BoxGeneration2, BoxGeneration3}
            setfield!(getfield(Wnew, x), :boxes, getfield(W, x).boxes)
            setfield!(getfield(Wnew, x), :probability, getfield(W, x).probability)
            setfield!(getfield(Wnew, x), :offset, getfield(W, x).offset)
            setfield!(getfield(Wnew, x), :dist, getfield(W, x).dist)
        elseif getfield(W, x) isa Random.AbstractRNG
            setfield!(Wnew, x, copy(getfield(W, x)))
        else
            # @warn "Got deep with $x::$(typeof(getfield(W, x))) in $(first(split(string(typeof(W)), '}')))"
            setfield!(Wnew, x, deepcopy(getfield(W, x)))
        end
    end
    # field u should be an alias for field W:
    if hasfield(typeof(W), :u)
        Wnew.u = Wnew.W
    end
    Wnew
end

function Base.copy(W::NoiseProcess)
    Wnew = NoiseProcess{isinplace(W)}(W.curt, W.curW, W.curZ, W.dist, W.bridge;
        rswm = W.rswm, save_everystep = W.save_everystep,
        rng = W.rng, covariance = W.covariance,
        reset = W.reset, reseed = W.reseed,
        continuous = W.continuous, cache = W.cache)
    copy!(Wnew, W)
end

function Base.copy(W::SimpleNoiseProcess)
    Wnew = SimpleNoiseProcess{isinplace(W)}(W.curt, W.curW, W.curZ, W.dist, W.bridge;
        save_everystep = W.save_everystep,
        rng = W.rng, covariance = W.covariance,
        reset = W.reset, reseed = W.reseed)
    copy!(Wnew, W)
end

function Base.copy(W::Union{NoiseWrapper, NoiseGrid, NoiseApproximation,
        VirtualBrownianTree, BoxWedgeTail})
    Wnew = typeof(W)((getfield(W, x) for x in fieldnames(typeof(W)))...)
    copy!(Wnew, W)
end

function Base.copy(W::NoiseFunction)
    Wnew = NoiseFunction{isinplace(W)}(W.t0, W.W, W.Z; noise_prototype = W.curW,
        reset = W.reset)
    copy!(Wnew, W)
end

function Base.copy(W::NoiseTransport)
    Wnew = NoiseTransport{isinplace(W)}(W.t0, W.W, W.RV, W.rv, W.Z;
        rng = W.rng,
        reset = W.reset, reseed = W.reseed,
        noise_prototype = W.curW)
    copy!(Wnew, W)
end
