function save_noise!(W::BoxWedgeTail)
    if W.t != W.curt
        push!(W.W, copy(W.curW))
        push!(W.t, copy(W.curt))
        push!(W.A, copy(W.curA))
        if W.Z != nothing
            push!(W.Z, copy(W.curZ))
        end
    end
    return nothing
end

@inline function accept_step!(W::BoxWedgeTail, dt, u, p, setup_next = true)

    W.curt += W.dt
    W.iter += 1

    if isinplace(W)
        @.. W.curW += W.dW
        if W.Z != nothing
            @.. W.curZ += W.dZ
        end
    else
        W.curW += W.dW
        if W.Z != nothing
            W.curZ += W.dZ
        end
    end

    if W.save_everystep
        push!(W.W, copy(W.curW))
        push!(W.t, copy(W.curt))
        push!(W.A, copy(W.curA))
        if W.Z != nothing
            push!(W.Z, copy(W.curZ))
        end
    end

    W.dt = dt #dtpropose
    # Setup next step
    if setup_next
        setup_next_step!(W::BoxWedgeTail, u, p)
    end
    return nothing
end

@inline function setup_next_step!(W::BoxWedgeTail, u, p)
    calculate_step!(W, W.dt, u, p)
    return nothing
end

@inline function calculate_step!(W::BoxWedgeTail, dt, u, p)
    if isinplace(W)
        sample_distribution(W.dW, W, dt, u, p, W.curt, W.rng)
        if W.Z != nothing
            W.dist(W.dZ, W, dt, u, p, W.curt, W.rng)
        end
    else
        W.dW, W.curA = sample_distribution(W, dt, u, p, W.curt, W.rng)
        if W.Z != nothing
            W.dZ = W.dist(W, dt, u, p, W.curt, W.rng)
        end
    end
    W.dt = dt
    return nothing
end

@inline function reject_step!(W::BoxWedgeTail, dtnew, u, p)
    error("BoxWedgeTail cannot be used with adaptivity rejections")
end

@inline function interpolate!(W::BoxWedgeTail, u, p, t)
    if sign(W.dt) * t > sign(W.dt) * W.t[end] # Steps past W
        dt = t - W.t[end]
        if isinplace(W)
            sample_distribution(W.dW, W, dt, u, p, t, W.rng)
            W.curW .+= W.dW
            if W.Z != nothing
                W.dist(W.dZ, W, dt, u, p, t, W.rng)
                W.curZ .+= W.dZ
            end
        else
            W.dW, W.curA = sample_distribution(W, dt, u, p, t, W.rng)
            W.curW += W.dW
            if W.Z != nothing
                W.dZ = W.dist(W, dt, u, p, t, W.rng)
                W.curZ += W.dZ
            end
        end
        out1 = copy(W.curW)
        if W.save_everystep
            push!(W.t, t)
            push!(W.W, out1)
        end
        if W.Z != nothing
            out2 = copy(W.curZ)
            if W.save_everystep
                push!(W.Z, out2)
            end
        else
            out2 = nothing
        end
        return out1, out2
    else # Bridge
        i = searchsortedfirst(W.t, t)
        if t == W.t[i]
            if isinplace(W)
                W.curW .= W.W[i]
            else
                W.curW = W.W[i]
            end
            if W.Z != nothing
                if isinplace(W)
                    W.curZ .= W.Z[i]
                else
                    W.curZ = W.Z[i]
                end
                return copy(W.curW), copy(W.curZ)
            else
                return copy(W.curW), nothing
            end
        else
            W0, Wh = W.W[i-1], W.W[i]
            if W.Z != nothing
                Z0, Zh = W.Z[i-1], W.Z[i]
            end
            h = W.t[i] - W.t[i-1]
            q = (t - W.t[i-1]) / h
            if isinplace(W)
                new_curW = similar(W.dW)
                W.bridge(new_curW, W, W0, Wh, q, h, u, p, t, W.rng)
                if iscontinuous(W)
                    @. new_curW += (1 - q) * W0
                else
                    @. new_curW += W0
                end
                if W.Z != nothing
                    new_curZ = similar(W.dZ)
                    W.bridge(new_curZ, W, Z0, Zh, q, h, u, p, t, W.rng)
                    if iscontinuous(W)
                        @. new_curZ += (1 - q) * Z0
                    else
                        @. new_curZ += Z0
                    end
                else
                    new_curZ = nothing
                end
            else
                new_curW = W.bridge(W, W0, Wh, q, h, u, p, t, W.rng)
                if iscontinuous(W)
                    # This should actually be based on the function for computing the mean
                    # flow of the noise process, but for now we'll just handle Wiener and
                    # Poisson
                    new_curW += (1 - q) * W0
                else
                    new_curW += W0
                end
                if W.Z != nothing
                    new_curZ = W.bridge(W, Z0, Zh, q, h, u, p, t, W.rng)
                    if iscontinuous(W)
                        new_curZ += (1 - q) * Z0
                    else
                        new_curZ += Z0
                    end
                else
                    new_curZ = nothing
                end
            end
            W.curW = new_curW
            if W.save_everystep
                insert!(W.W, i, new_curW)
                insert!(W.t, i, t)
            end
            if W.Z != nothing
                W.curZ = new_curZ
                if W.save_everystep
                    insert!(W.Z, i, new_curZ)
                end
            end
            return new_curW, new_curZ
        end
    end
end

@inline function interpolate!(out1, out2, W::BoxWedgeTail, u, p, t)
    if sign(W.dt) * t > sign(W.dt) * W.t[end] # Steps past W
        dt = t - W.t[end]
        sample_distribution(W.dW, W, dt, u, p, t, W.rng)
        out1 .+= W.dW
        if W.Z != nothing
            W.dist(W.dZ, W, dt, u, p, t, W.rng)
            out2 .+= W.dZ
        end
        if W.save_everystep
            push!(W.t, t)
            push!(W.W, copy(out1))
            if W.Z != nothing
                push!(W.Z, copy(out2))
            end
        end
    else # Bridge
        i = searchsortedfirst(W.t, t)
        if t == W.t[i]
            out1 .= W.W[i]
            if W.Z != nothing
                out2 .= W.Z[i]
            end
        else
            W0, Wh = W.W[i-1], W.W[i]
            if W.Z != nothing
                Z0, Zh = W.Z[i-1], W.Z[i]
            end
            h = W.t[i] - W.t[i-1]
            q = (t - W.t[i-1]) / h
            W.bridge(out1, W, W0, Wh, q, h, u, p, t, W.rng)
            out1 .+= (1 - q) * W0
            if W.Z != nothing
                W.bridge(out2, W, Z0, Zh, q, h, u, p, t, W.rng)
                out2 .+= (1 - q) * Z0
            end
            W.curW .= out1
            if W.save_everystep
                insert!(W.W, i, copy(out1))
                insert!(W.t, i, t)
            end
            if W.Z != nothing
                W.curZ .= out2
                if W.save_everystep
                    insert!(W.Z, i, copy(out2))
                end
            end
        end
    end
end

# joint probability densitity function
function joint_density_function(r, a, rtol)
    integral, err = QuadGK.quadgk(
        x -> r / π * (x / sinh(x)) * exp(-r^2 * x / (2 * tanh(x))) * cos(a * x),
        0,
        Inf,
        rtol = rtol,
    )
    return integral
end

# boxes
abstract type AbstractBoxGeneration end

mutable struct BoxGeneration1{boxesType,pType,offType,distType} <: AbstractBoxGeneration
    boxes::boxesType
    probability::pType
    offset::offType
    dist::distType
end

struct BoxGeneration2{boxesType,pType,offType,distType} <: AbstractBoxGeneration
    boxes::boxesType
    probability::pType
    offset::offType
    dist::distType
end

struct BoxGeneration3{boxesType,pType,offType,distType} <: AbstractBoxGeneration
    boxes::boxesType
    probability::pType
    offset::offType
    dist::distType
end

function generate_boxes1(densf, Δr, Δa, Δz, rM, aM, offset = nothing, scale = 1)
    boxes = Array{typeof(Δr),1}[]
    probability = Vector{typeof(Δr)}(undef, 0)
    rs = collect(zero(rM):Δr:rM)
    as = collect(zero(aM):Δa:aM)

    if offset === nothing
        offset = zeros(typeof(Δr), Int(rM / Δr), Int(aM / Δa))
    end

    for (i, r) in enumerate(rs[1:end-1])
        for (j, a) in enumerate(as[1:end-1])
            # height of the box
            z = Δz
            if offset !== nothing
                indx1 = (i - 1) * scale + one(scale)
                indx2 = (j - 1) * scale + one(scale)
                z += offset[indx1, indx2]
            end

            counter = 0

            # check for complete inclusion of box in volume
            while (
                z <= densf(r, a) &&
                z <= densf(r + Δr, a) &&
                z <= densf(r, a + Δa) &&
                z <= densf(r + Δr, a + Δa)
            )
                counter += one(counter)
                z += Δz
            end
            z -= Δz
            # store position of top corner of box and width
            if counter != 0
                push!(boxes, [r, a, Δr, Δa])
                push!(probability, counter * 2 * Δr * Δa * Δz)
            end

            if offset !== nothing
                offset[indx1:indx1+(scale-one(scale)), indx2:indx2+(scale-one(scale))] .= z
            end

        end
    end
    return boxes, probability, offset
end

function generate_boxes2(
    densf,
    Δrmin,
    Δamin,
    Δzmin,
    Δrmax,
    Δamax,
    Δzmax,
    rM,
    aM,
    offset = nothing,
)
    boxes = Array{typeof(Δrmin),1}[]
    probability = Vector{typeof(Δrmin)}(undef, 0)
    # start with largest possible size, then subsequently decrease size and fill remaining space
    Δr = Δrmax
    Δa = Δamax
    Δz = Δzmax

    if offset === nothing
        offset = zeros(typeof(Δrmin), Int(rM / Δrmin), Int(aM / Δamin))
    end

    while (Δr >= Δrmin && Δa >= Δamin && Δz >= Δzmin)
        scale = Int(Δr / Δrmin)
        boxestmp, probabilitytmp = generate_boxes1(densf, Δr, Δa, Δz, rM, aM, offset, scale)
        boxes = cat(boxes, boxestmp, dims = 1)
        probability = cat(probability, probabilitytmp, dims = 1)
        Δr = Δr / 2
        Δa = Δa / 2
        Δz = Δz / 2
    end
    # merge boxes
    return boxes, probability, offset
end

function generate_boxes3(densf, Δr, Δa, Δz, rM, aM, offset = nothing, scale = 1)
    boxes = Array{typeof(Δr),1}[]
    probability = Vector{typeof(Δr)}(undef, 0)
    rs = collect(zero(rM):Δr:rM)
    as = collect(zero(aM):Δa:aM)

    if offset === nothing
        offset = zeros(typeof(Δr), Int(rM / Δr), Int(aM / Δa))
    end

    for (i, r) in enumerate(rs[1:end-1])
        for (j, a) in enumerate(as[1:end-1])
            # height of the box
            z = Δz
            if offset !== nothing
                indx1 = (i - 1) * scale + one(scale)
                indx2 = (j - 1) * scale + one(scale)
                z += offset[indx1, indx2]
            end
            # check for complete inclusion of box in volume
            while (
                z <= densf(r, a) &&
                z <= densf(r + Δr, a) &&
                z <= densf(r, a + Δa) &&
                z <= densf(r + Δr, a + Δa)
            )
                # store position of top corner of box and width
                push!(boxes, [r, a, Δr, Δa])
                # probability = 2*Δr*Δa*Δz for all small boxes, push!(probability, 1)
                push!(probability, 2 * Δr * Δa * Δz)

                z += Δz
            end
            z -= Δz
            if offset !== nothing
                offset[indx1:indx1+(scale-one(scale)), indx2:indx2+(scale-one(scale))] .= z
            end
        end
    end
    return boxes, probability, offset
end


# entropy of partition
function entropy(probs)
    TotalVolume = sum(probs)
    qi = probs ./ TotalVolume
    return -sum(qi .* log2.(qi))
end


# sampling from boxes
function sample_box(W::BoxWedgeTail, Boxes::AbstractBoxGeneration)
    indx = rand(W.rng, Boxes.dist)
    # boxes store r, a, Δr, Δa 
    ri, ai, Δr, Δa = Boxes.boxes[indx]
    DU = Distributions.Product(Distributions.Uniform.([ri, ai], [ri + Δr, ai + Δa]))

    r, a = rand(W.rng, DU)
    return r, a
end


# Wedges
struct Wedges{boxesType,pType,offType,distType} <: AbstractBoxGeneration
    boxes::boxesType
    probability::pType
    fvalues::offType
    dist::distType
end

# linear approximation for squeezing method
function linear_interpolation_wedges(fij, fij2, fij3, fij4, r, a, ri, ai, Δr, Δa)
    t = (r - ri) / Δr
    u = (a - ai) / Δa
    return (one(t) - t) * (one(u) - u) * fij +
           t * (one(u) - u) * fij2 +
           (one(t) - t) * u * fij3 +
           t * u * fij4
end

function contrained_optimization_problem(densf, fij, fij2, fij3, fij4, ri, ai, Δr, Δa)
    function difference(x)
        densf(x[1], x[2]) -
        linear_interpolation_wedges(fij, fij2, fij3, fij4, x[1], x[2], ri, ai, Δr, Δa)
    end
    ϵijmax = Optim.optimize(
        x -> -difference(x),
        [ri, ai],
        [ri + Δr, ai + Δa],
        [ri + Δr / 2, ai + Δa / 2],
        Optim.Fminbox(Optim.NelderMead()),
    )
    ϵijmin = Optim.optimize(
        x -> difference(x),
        [ri, ai],
        [ri + Δr, ai + Δa],
        [ri + Δr / 2, ai + Δa / 2],
        Optim.Fminbox(Optim.NelderMead()),
    )
    return ϵijmin.minimum, ϵijmax.minimum
end

function generate_wedges(densf, Δr, Δa, Δz, rM, aM, offset, sqeezing)
    boxes = Array{typeof(Δr),1}[]
    probability = Vector{typeof(Δr)}(undef, 0)
    f = Array{typeof(Δr),1}[]
    rs = collect(zero(rM):Δr:rM)
    as = collect(zero(aM):Δa:aM)
    for (i, r) in enumerate(rs[1:end-1])
        for (j, a) in enumerate(as[1:end-1])
            # base height hij of the box enclosing the wedge
            hij = offset[i, j]
            # max height of the box f̃_ij
            fij = densf(r, a) # f_{ij}
            fij2 = densf(r + Δr, a) # f_{i+1,j}
            fij3 = densf(r, a + Δa) # f_{i,j+1}
            fij4 = densf(r + Δr, a + Δa) # f_{i+1,j+1}

            f̃ij = max(fij, fij2, fij3, fij4)

            # store position of top corner of box and width
            if sqeezing
                ϵijmin, ϵijmax = contrained_optimization_problem(
                    densf,
                    fij,
                    fij2,
                    fij3,
                    fij4,
                    r,
                    a,
                    Δr,
                    Δa,
                )
                push!(boxes, [f̃ij, hij, abs(ϵijmin), abs(ϵijmax), r, a])
                # store PDF values
                push!(f, [fij, fij2, fij3, fij4])
            else
                push!(boxes, [f̃ij, hij, zero(hij), zero(hij), r, a])
            end

            # probability = 2*Δr*Δa*Δz
            push!(probability, 2 * Δr * Δa * (f̃ij - hij))
        end
    end
    return boxes, probability, f
end


# sampling from wedges
function sample_wedge(W::BoxWedgeTail, wedges::Wedges)
    indx = rand(W.rng, wedges.dist)
    # wedges store f̃ij, hij, ϵijmin, ϵijmax, r, a, Δr
    f̃ij, hij, ϵijmin, ϵijmax, ri, ai = wedges.boxes[indx]
    DU = Distributions.Product(
        Distributions.Uniform.([ri, ai, hij], [ri + W.Δr, ai + W.Δa, f̃ij]),
    )
    if W.sqeezing
        fij, fij2, fij3, fij4 = wedges.fvalues[indx]
        while true
            r, a, z = rand(W.rng, DU)
            flin =
                linear_interpolation_wedges(fij, fij2, fij3, fij4, r, a, ri, ai, W.Δr, W.Δa)
            if z > flin + ϵijmax
                continue
            elseif z < flin - ϵijmin
                break
            end

            if z <= W.jpdf(r, a)
                break
            end
        end
    else
        r, a, z = rand(W.rng, DU)
        while z > W.jpdf(r, a)
            r, a, z = rand(W.rng, DU)
        end
    end
    return r, a
end


# Tail approximation
abstract type AbstractTail end

struct Tail1{pType,distType,pdfType,cType} <: AbstractTail
    p::pType
    dist::distType
    candpdf::pdfType
    cvalue::cType # success probability 1/c.

    function Tail1(rM, aM)
        p = convert(typeof(rM), 0.0002982405821953734)
        candpdf = (r, a) -> exp(-r^2 / 2)
        dist1 = Distributions.truncated(
            Distributions.Normal(zero(rM), one(rM)),
            rM,
            12 * one(rM),
        )
        dist2 = Distributions.Uniform(zero(aM), aM)

        dist = Distributions.Product([dist1, dist2])

        c = convert(typeof(rM), 2.3)

        new{typeof(p),typeof(dist),typeof(candpdf),typeof(c)}(p, dist, candpdf, c)
    end
end

struct Tail2{pType,distType,pdfType,cType,boundsType} <: AbstractTail
    p::pType
    dist::distType
    candpdf::pdfType
    cvalue::cType # success probability 1/c.
    alower::boundsType
    aupper::boundsType

    function Tail2(rM, aM)
        p = convert(typeof(rM), 3.648002197730662e-5)
        candpdf = (r, a) -> 3 // 2 * exp(-r^2 / 2 - 3 * a^2 / (2 * r^2))
        # store only distribution of r, when sampling re-create distribution for a|r
        dist = Distributions.truncated(
            Distributions.Normal(zero(rM), one(rM)),
            rM,
            8 * one(rM),
        )

        c = convert(typeof(rM), 2.0)

        alower = aM
        aupper = 2 * aM

        new{typeof(p),typeof(dist),typeof(candpdf),typeof(c),typeof(alower)}(
            p,
            dist,
            candpdf,
            c,
            alower,
            aupper,
        )
    end
end

struct Tail3{pType,distType,pdfType,cType} <: AbstractTail
    p::pType
    dist::distType
    candpdf::pdfType
    cvalue::cType # success probability 1/c.

    function Tail3(rM, aM)
        p = convert(typeof(rM), 0.0018026889941638695)
        candpdf = (r, a) -> exp(convert(typeof(rM), -pi / 2) * a)
        dist1 = Distributions.Uniform(2 * one(rM), rM)
        dist2 = Distributions.truncated(
            Distributions.Exponential(convert(typeof(rM), 2 / pi)),
            aM,
            8 * one(aM),
        )

        dist = Distributions.Product([dist1, dist2])

        c = convert(typeof(rM), 2.6)

        new{typeof(p),typeof(dist),typeof(candpdf),typeof(c)}(p, dist, candpdf, c)
    end
end

struct Tail4{pType,distType,pdfType,cType} <: AbstractTail
    p::pType
    dist::distType
    candpdf::pdfType
    cvalue::cType # success probability 1/c.

    function Tail4(rM, aM)
        p = convert(typeof(rM), 1.6145323459108908e-6)
        candpdf = (r, a) -> 15 * r * exp(convert(typeof(rM), -pi) * a)
        dist1 = Distributions.TriangularDist(zero(rM), one(rM) / 2, one(rM) / 2)
        dist2 = Distributions.truncated(
            Distributions.Exponential(convert(typeof(rM), 1 / pi)),
            aM,
            6 * one(aM),
        )

        dist = Distributions.Product([dist1, dist2])

        c = convert(typeof(rM), 2.6)

        new{typeof(p),typeof(dist),typeof(candpdf),typeof(c)}(p, dist, candpdf, c)
    end
end

struct Tail5{pType,distType,pdfType,cType} <: AbstractTail
    p::pType
    dist::distType
    candpdf::pdfType
    cvalue::cType # success probability 1/c.

    function Tail5(rM, aM)
        p = convert(typeof(rM), 1.9637270929071978e-5)
        candpdf = (r, a) -> 15 * r * exp(convert(typeof(rM), -2.8) * a)
        dist1 = Distributions.TriangularDist(one(rM) / 2, one(rM), one(rM))
        dist2 = Distributions.truncated(
            Distributions.Exponential(convert(typeof(rM), 1 / 2.8)),
            aM,
            6 * one(aM),
        )

        dist = Distributions.Product([dist1, dist2])

        c = convert(typeof(rM), 2.8)

        new{typeof(p),typeof(dist),typeof(candpdf),typeof(c)}(p, dist, candpdf, c)
    end
end

struct Tail6{pType,distType,pdfType,cType} <: AbstractTail
    p::pType
    dist::distType
    candpdf::pdfType
    cvalue::cType # success probability 1/c.

    function Tail6(rM, aM)
        p = convert(typeof(rM), 0.00012022274714710923)
        candpdf = (r, a) -> 25 * r * exp(convert(typeof(rM), -2.6) * a)
        dist1 = Distributions.TriangularDist(one(rM), 3 / 2 * one(rM), 3 / 2 * one(rM))
        dist2 = Distributions.truncated(
            Distributions.Exponential(convert(typeof(rM), 1 / 2.6)),
            aM,
            6 * one(aM),
        )

        dist = Distributions.Product([dist1, dist2])

        c = convert(typeof(rM), 3.0)

        new{typeof(p),typeof(dist),typeof(candpdf),typeof(c)}(p, dist, candpdf, c)
    end
end

struct Tail7{pType,distType,pdfType,cType} <: AbstractTail
    p::pType
    dist::distType
    candpdf::pdfType
    cvalue::cType # success probability 1/c.

    function Tail7(rM, aM)
        p = convert(typeof(rM), 0.00038595770876898326)
        candpdf = (r, a) -> 25 * r * exp(convert(typeof(rM), -2.4) * a)
        dist1 = Distributions.TriangularDist(3 / 2 * one(rM), 2 * one(rM), 2 * one(rM))
        dist2 = Distributions.truncated(
            Distributions.Exponential(convert(typeof(rM), 1 / 2.4)),
            aM,
            6 * one(aM),
        )

        dist = Distributions.Product([dist1, dist2])

        c = convert(typeof(rM), 3.2)

        new{typeof(p),typeof(dist),typeof(candpdf),typeof(c)}(p, dist, candpdf, c)
    end
end

struct Tail8{pType,distType,pdfType,cType} <: AbstractTail
    p::pType
    dist::distType
    candpdf::pdfType
    cvalue::cType # success probability 1/c.

    function Tail8(rM, aM)
        p = convert(typeof(rM), 6.566973361922312e-6)
        candpdf = (r, a) -> 40 * r * exp(convert(typeof(rM), -2.4) * a)
        dist1 = Distributions.TriangularDist(one(rM), 2 * one(rM), 2 * one(rM))
        dist2 = Distributions.truncated(
            Distributions.Exponential(convert(typeof(rM), 1 / 2.4)),
            6 * one(aM),
            8 * one(aM),
        )

        dist = Distributions.Product([dist1, dist2])

        c = convert(typeof(rM), 4.2)

        new{typeof(p),typeof(dist),typeof(candpdf),typeof(c)}(p, dist, candpdf, c)
    end
end

struct Tail9{pType,distType,pdfType,cType} <: AbstractTail
    p::pType
    dist::distType
    candpdf::pdfType
    cvalue::cType # success probability 1/c.

    function Tail9(rM, aM)
        p = convert(typeof(rM), 4.149989421681627e-6)
        candpdf = (r, a) -> 7 // 10 * exp(convert(typeof(rM), -pi / 2) * a)
        dist1 = Distributions.Uniform(2 * one(rM), 5 * one(rM))
        dist2 = Distributions.truncated(
            Distributions.Exponential(convert(typeof(rM), 2 / pi)),
            8 * one(aM),
            10 * one(aM),
        )

        dist = Distributions.Product([dist1, dist2])

        c = convert(typeof(rM), 2.4)

        new{typeof(p),typeof(dist),typeof(candpdf),typeof(c)}(p, dist, candpdf, c)
    end
end

struct TailApproxs{
    T1<:Tail1,
    T2<:Tail2,
    T3<:Tail3,
    T4<:Tail4,
    T5<:Tail5,
    T6<:Tail6,
    T7<:Tail7,
    T8<:Tail8,
    T9<:Tail9,
    distType,
}
    tail1::T1
    tail2::T2
    tail3::T3
    tail4::T4
    tail5::T5
    tail6::T6
    tail7::T7
    tail8::T8
    tail9::T9
    dist::distType

    function TailApproxs(rM, aM)
        tail1 = Tail1(rM, aM)
        tail2 = Tail2(rM, aM)
        tail3 = Tail3(rM, aM)
        tail4 = Tail4(rM, aM)
        tail5 = Tail5(rM, aM)
        tail6 = Tail6(rM, aM)
        tail7 = Tail7(rM, aM)
        tail8 = Tail8(rM, aM)
        tail9 = Tail9(rM, aM)

        probability = [
            tail1.p,
            tail2.p,
            tail3.p,
            tail4.p,
            tail5.p,
            tail6.p,
            tail7.p,
            tail8.p,
            tail9.p,
        ]

        dist = Distributions.Categorical(probability / sum(probability))

        new{
            typeof(tail1),
            typeof(tail2),
            typeof(tail3),
            typeof(tail4),
            typeof(tail5),
            typeof(tail6),
            typeof(tail7),
            typeof(tail8),
            typeof(tail9),
            typeof(dist),
        }(
            tail1,
            tail2,
            tail3,
            tail4,
            tail5,
            tail6,
            tail7,
            tail8,
            tail9,
            dist,
        )
    end
end

# sampling from tail approximations by rejection sampling
function sample_tail(W::BoxWedgeTail, tails::TailApproxs)
    indx = rand(W.rng, tails.dist)
    # decide from which tail region to sample
    if indx == 1
        T = tails.tail1
    elseif indx == 2
        T = tails.tail2
    elseif indx == 3
        T = tails.tail3
    elseif indx == 4
        T = tails.tail4
    elseif indx == 5
        T = tails.tail5
    elseif indx == 6
        T = tails.tail6
    elseif indx == 7
        T = tails.tail7
    elseif indx == 8
        T = tails.tail8
    else
        T = tails.tail9
    end

    r, a = sample_tail(W.rng, W.jpdf, T)
end

function sample_tail(rng, densf, T::Tail2)

    local r, a

    while true
        # draw uniform
        U = rand(rng)
        # draw a candidate X∼g from the candidate density
        # draw first the r value
        r = rand(rng, T.dist)
        # generate a from the conditional distribution
        dista = Distributions.truncated(
            Distributions.Normal(zero(r), sqrt(r / 3)),
            T.alower,
            T.aupper,
        )
        a = rand(rng, dista)
        if U <= densf(r, a) / (T.cvalue * T.candpdf(r, a))
            break
        end
    end

    return r, a
end

function sample_tail(rng, densf, T)
    local r, a
    while true
        # draw uniform
        U = rand(rng)
        # draw a candidate X∼g from the candidate density
        r, a = rand(rng, T.dist)
        if U <= densf(r, a) / (T.cvalue * T.candpdf(r, a))
            break
        end
    end

    return r, a
end


function sample_distribution(dW, W::BoxWedgeTail, dt, u, p, t, rng)
    # decide if sample should come from boxes, wedges, or tails
    indx = rand(rng, W.distBWT)

    # draw (r,a) pair
    if indx == 1
        r, a = DiffEqNoiseProcess.sample_box(W, W.boxes)
    elseif indx == 2
        r, a = DiffEqNoiseProcess.sample_wedge(W, W.wedges)
    else
        r, a = DiffEqNoiseProcess.sample_tail(W, W.tails)
    end

    # get dWs from r and scale by sqrt(dt)
    θ = rand(rng, W.distΠ)
    sqabsdt = @fastmath sqrt(abs(dt))
    @. dW = sqabsdt * r * [cos(θ), sin(θ)]

    # random sign for A and scale by dt
    W.curA = rand(rng, (-1, 1)) * a * dt

    return nothing
end

function sample_distribution(W::BoxWedgeTail, dt, u, p, t, rng)
    # decide if sample should come from boxes, wedges, or tails
    indx = rand(rng, W.distBWT)

    # draw (r,a) pair
    if indx == 1
        r, a = DiffEqNoiseProcess.sample_box(W, W.boxes)
    elseif indx == 2
        r, a = DiffEqNoiseProcess.sample_wedge(W, W.wedges)
    else
        r, a = DiffEqNoiseProcess.sample_tail(W, W.tails)
    end

    # get dWs from r and scale by sqrt(dt)
    θ = rand(rng, W.distΠ)
    sqabsdt = @fastmath sqrt(abs(dt))
    dW = sqabsdt * r * [cos(θ), sin(θ)]

    # random sign for A and scale by dt
    A = rand(rng, (-1, 1)) * a * dt

    return dW, A
end
