function save_noise!(W::VirtualBrownianTree)
end

function interpolate!(W::VirtualBrownianTree, u, p, t)
    ts, timeseries, timeseries2 = W.t, W.W, W.Z

    if t > ts[end] ≥ ts[1] || t < ts[end] < ts[1]
        error("Bridge cannot extrapolate past the final timepoint. Build a longer VBT to cover the integration.")
    end

    if t < ts[1] ≤ ts[end] || t > ts[1] ≥ ts[end]
        error("Bridge cannot extrapolate before the first timepoint. Build a longer VBT to cover the integration.")
    end

    tdir = sign(ts[end] - ts[1])

    if t isa Union{Rational, Integer}
        @inbounds i = searchsortedfirst(ts, t, rev = tdir < 0) # It's in the interval ts[i-1] to ts[i]
    else
        @inbounds i = searchsortedfirst(ts, t - tdir * 10eps(typeof(t)), rev = tdir < 0)
    end

    @inbounds if (t isa Union{Rational, Integer} && ts[i] == t) ||
                 (isapprox(t, ts[i]; atol = 100eps(typeof(t)), rtol = 100eps(t)))
        val1 = timeseries[i]
        W.dZ !== nothing ? val2 = timeseries2[i] : val2 = nothing
    elseif ts[i - 1] == t # Can happen if it's the first value!
        val1 = timeseries[i - 1]
        W.dZ !== nothing ? val2 = timeseries2[i - 1] : val2 = nothing
    else
        if W.dZ !== nothing
            val1, val2 = search_VBT(t, W.seeds[i - 1], ts[i - 1], ts[i], timeseries[i - 1],
                                    timeseries[i],
                                    timeseries2[i - 1], timeseries2[i], W, W.rng)
        else
            val1, val2 = search_VBT(t, W.seeds[i - 1], ts[i - 1], ts[i], timeseries[i - 1],
                                    timeseries[i],
                                    nothing, nothing, W, W.rng)
        end
    end
    val1, val2
end

function interpolate!(out1, out2, W::VirtualBrownianTree, u, p, t)
    ts, timeseries, timeseries2 = W.t, W.W, W.Z
    sign(W.dt) * t > sign(W.dt) * ts[end] &&
        error("Bridge cannot extrapolate past the final timepoint. Build a longer VBT to cover the integration.")
    sign(W.dt) * t < sign(W.dt) * ts[1] &&
        error("Bridge cannot extrapolate before the first timepoint. Build a longer VBT to cover the integration.")
    tdir = sign(ts[end] - ts[1])

    if t isa Union{Rational, Integer}
        @inbounds i = searchsortedfirst(ts, t, rev = tdir < 0) # It's in the interval ts[i-1] to ts[i]
    else
        @inbounds i = searchsortedfirst(ts, t - tdir * 10eps(typeof(t)), rev = tdir < 0)
    end

    @inbounds if (t isa Union{Rational, Integer} && ts[i] == t) ||
                 (isapprox(t, ts[i]; atol = 100eps(typeof(t)), rtol = 100eps(t)))
        copyto!(out1, timeseries[i])
        W.dZ !== nothing && copyto!(out2, timeseries2[i])
    elseif ts[i - 1] == t # Can happen if it's the first value!
        copyto!(out1, timeseries[i - 1])
        W.dZ !== nothing && copyto!(out2, timeseries2[i - 1])
    else
        if W.dZ !== nothing
            search_VBT!(out1, out2, t, seeds[i - 1], ts[i - 1], ts[i], timeseries[i - 1],
                        timeseries[i],
                        timeseries2[i - 1], timeseries2[i], W, W.rng)
        else
            search_VBT!(out1, out2, t, seeds[i - 1], ts[i - 1], ts[i], timeseries[i - 1],
                        timeseries[i],
                        nothing, nothing, W, W.rng)
        end
    end
    return nothing
end

function calculate_step!(W::VirtualBrownianTree, dt, u, p)
    t = W.curt + dt
    if typeof(t) <: AbstractFloat && abs(t - W.t[end]) < 100eps(typeof(dt))
        t = W.t[end]
    end
    if isinplace(W)
        interpolate!(W.dW, W.dZ, W, u, p, t)
        W.dW .-= W.curW
        if W.dZ !== nothing
            W.dZ .-= W.curZ
        end
    else
        new_W, new_Z = interpolate!(W, u, p, t)
        W.dW = new_W - W.curW
        if W.dZ !== nothing
            W.dZ = new_Z - W.curZ
        end
    end
    W.dt = dt
    return nothing
end

function accept_step!(W::VirtualBrownianTree, dt, u, p, setup_next = true)
    W.step_setup == false && error("Stepped past the defined domain for the VBT.")

    if isinplace(W)
        W.curW .+= W.dW
    else
        W.curW += W.dW
    end
    W.curt += W.dt
    if W.dZ !== nothing
        if isinplace(W)
            W.curZ .+= W.dZ
        else
            W.curZ += W.dZ
        end
    end

    W.dt = dt #dtpropose
    t_min, t_max = extrema((W.t[begin], W.t[end]))
    t = W.curt + W.dt
    if (W.dt isa Union{Rational, Integer})
        if t < t_min || t > t_max
            setup_next = false
            W.step_setup = false
        end
    else
        if t < t_min - 10eps(typeof(dt)) || t > t_max + 10eps(typeof(dt))
            setup_next = false
            W.step_setup = false
        end
    end

    if setup_next
        calculate_step!(W, dt, u, p)
    end
    return nothing
end

function reject_step!(W::VirtualBrownianTree, dtnew, u, p)
    calculate_step!(W, dtnew, u, p)
    return nothing
end

function setup_next_step!(W::VirtualBrownianTree, u, p)
    calculate_step!(W, W.dt, u, p)
    return nothing
end

# split seeds

# for counter-based RNGs
function split_VBT_seed(rng::Random123.AbstractR123, parent_seed, current_depth, Nt)

    # seed left
    seed_l = convert(typeof(parent_seed), parent_seed - (Nt - 1) / 2^(current_depth + 1))
    # seed right
    seed_r = convert(typeof(parent_seed), parent_seed + (Nt - 1) / 2^(current_depth + 1))

    RandomNumbers.Random123.set_counter!(rng, parent_seed)
    return seed_l, seed_r, parent_seed
end

# create the cache of the VBT with depth tree_depth
function create_VBT_cache(bridge, t0, W0, Z0, tend, Wend, Zend, rng::Random123.AbstractR123,
                          tree_depth, search_depth)
    # total number of cached time steps and W values
    Nt = Int(2^search_depth + 1)

    ts = [t0, tend]
    Ws = [W0, Wend]

    dW = Wend - W0 # for input of bridges

    if Z0 !== nothing
        Zs = [Z0, Zend]
    else
        Zs = [nothing]
    end

    seeds = [convert(typeof(rng.ctr1), rng.ctr1 + (Nt + 1) / 2)]

    q = 1 // 2

    for level in 1:Int(tree_depth)
        new_ts = Vector{typeof(t0)}(undef, 0)
        new_Ws = Vector{typeof(W0)}(undef, 0)
        if Z0 !== nothing
            new_Zs = Vector{typeof(Z0)}(undef, 0)
        else
            new_Zs = [nothing]
        end
        new_seeds = Vector{typeof(rng.ctr1)}(undef, 0)
        for (i, parent) in enumerate(seeds)
            seed_l, seed_r, seed_v = split_VBT_seed(rng, parent, level, Nt)
            append!(new_seeds, [seed_l, seed_r])

            t0, t1 = ts[i], ts[i + 1]
            W0tmp, W1tmp = Ws[i], Ws[i + 1]
            t = (t0 + t1) / 2

            h = t1 - t0
            #q = (t-t0)/h # == 1//2 defined above

            w = bridge(dW, nothing, W0tmp, W1tmp, q, h, nothing, nothing, nothing, rng)
            if Z0 !== nothing
                Z0tmp, Z1tmp = Zs[i], Zs[i + 1]
                z = bridge(dW, nothing, Z0tmp, Z1tmp, q, h, nothing, nothing, nothing, rng)
                append!(new_Zs, Z0)
                append!(new_Zs, z)
            end

            append!(new_ts, [t0, t])
            append!(new_Ws, [W0tmp, w])
        end
        push!(new_ts, tend)
        push!(new_Ws, Wend)

        if Z0 !== nothing
            push!(new_Zs, Zend)
            Zs = new_Zs
        end

        ts = new_ts
        Ws = new_Ws

        seeds = new_seeds
    end

    return ts, Ws, Zs, seeds
end

function search_VBT(t, seed, t0, t1, W0, W1, Z0, Z1, W::VirtualBrownianTree,
                    rng::Random123.AbstractR123)
    Nt = Int(2^W.search_depth + 1)
    depth = Int(W.tree_depth + 1)
    seed_l, seed_r, seed_v = split_VBT_seed(rng, seed, depth, Nt)

    tmid = (t0 + t1) / 2
    h = t1 - t0
    q = 1 // 2

    out1 = W.bridge(W.dW, nothing, W0, W1, q, h, nothing, nothing, nothing, rng)

    if Z0 !== nothing
        out2 = W.bridge(W.dW, nothing, Z0, Z1, q, h, nothing, nothing, nothing, rng)
    else
        out2 = nothing
    end

    # create a buffer
    W0tmp, W1tmp = copy(W0), copy(W1)
    if Z0 !== nothing
        Z0tmp, Z1tmp = copy(Z0), copy(Z1)
    end

    while (abs(t - tmid) > W.atol && depth < W.search_depth)
        depth += 1
        if t < tmid
            t1 = tmid
            W0tmp, W1tmp = W0tmp, out1
            seed_v = seed_l
            if Z0 !== nothing
                Z0tmp, Z1tmp = Z0tmp, out2
            end
        else
            t0 = tmid
            W0tmp, W1tmp = out1, W1tmp
            seed_v = seed_r
            if Z0 !== nothing
                Z0tmp, Z1tmp = out2, Z1tmp
            end
        end

        seed_l, seed_r, seed_v = split_VBT_seed(rng, seed_v, depth, Nt)
        tmid = (t0 + t1) / 2
        h = t1 - t0

        out1 = W.bridge(W.dW, nothing, W0tmp, W1tmp, q, h, nothing, nothing, nothing, rng)
        if Z0 !== nothing
            out2 = W.bridge(W.dW, nothing, Z0tmp, Z1tmp, q, h, nothing, nothing, nothing,
                            rng)
        end
    end

    return out1, out2
end

function search_VBT!(out1, out2, t, seed, t0, t1, W0, W1, Z0, Z1, W::VirtualBrownianTree,
                     rng::Random123.AbstractR123)
    Nt = Int(2^W.search_depth + 1)
    depth = Int(W.tree_depth + 1)
    seed_l, seed_r, seed_v = split_VBT_seed(rng, seed, depth, Nt)

    tmid = (t0 + t1) / 2
    h = t1 - t0
    q = 1 // 2

    # create a buffer
    copyto!(W.W0tmp, W0)
    copyto!(W.W1tmp, W1)
    if Z0 !== nothing
        copyto!(W.Z0tmp, Z0)
        copyto!(W.Z1tmp, Z1)
    end

    W.bridge(out1, W.dW, nothing, W.W0tmp, W.W1tmp, q, h, nothing, nothing, nothing, rng)

    if Z0 !== nothing
        W.bridge(out2, W.dW, nothing, W.Z0tmp, W.Z1tmp, q, h, nothing, nothing, nothing,
                 rng)
    end

    while (abs(t - tmid) > W.atol && depth < W.search_depth)
        depth += 1
        if t < tmid
            t1 = tmid
            copyto!(W.W1tmp, out1)
            seed_v = seed_l
            if Z0 !== nothing
                copyto!(W.Z1tmp, out2)
            end
        else
            t0 = tmid
            copyto!(W.W0tmp, out1)
            seed_v = seed_r
            if Z0 !== nothing
                copyto!(W.Z0tmp, out2)
            end
        end

        seed_l, seed_r, seed_v = split_VBT_seed(rng, seed_v, depth, Nt)
        tmid = (t0 + t1) / 2
        h = t1 - t0

        W.bridge(out1, W.dW, nothing, W.W0tmp, W.W1tmp, q, h, nothing, nothing, nothing,
                 rng)
        if Z0 !== nothing
            W.bridge(out2, W.dW, nothing, W.Z0tmp, W.Z1tmp, q, h, nothing, nothing, nothing,
                     rng)
        end
    end

    return nothing
end
