@testset "RSwM rejection correctness" begin
    using DiffEqNoiseProcess, Test, Random, Statistics

    # A NoiseProcess driven through a FIXED accept/reject pattern must still produce
    # a genuine Brownian path: increments over the accepted grid are independent and
    # N(0, Δt). Because the pattern is fixed, the grid is identical across seeds, so
    # the increments can be compared statistically.
    #
    # Rejection factors are expressed relative to W.dt rather than as absolute step
    # sizes, since RSwM1's setup_next_step! overrides W.dt from the stack.

    # (:reject, q) shrinks to q*W.dt; (:accept, g) proposes g*W.dt for the next step
    function make_schedule(; nacc = 8, seed = 1, preject = 0.55)
        rng = Xoshiro(seed)
        ops = Tuple{Symbol, Float64}[]
        for _ in 1:nacc
            while rand(rng) < preject
                push!(ops, (:reject, 0.2 + 0.7 * rand(rng)))
            end
            push!(ops, (:accept, 1.0 + 1.5 * rand(rng)))
        end
        return ops
    end

    function build(alg, iip, n, withZ, seed)
        return if iip
            WienerProcess!(
                0.0, zeros(n), withZ ? zeros(n) : nothing;
                rswm = RSWM(adaptivealg = alg), rng = Xoshiro(seed)
            )
        else
            WienerProcess(
                0.0, 0.0, withZ ? 0.0 : nothing;
                rswm = RSWM(adaptivealg = alg), rng = Xoshiro(seed)
            )
        end
    end

    function run_path(ops, alg, iip, n, withZ, seed; h0 = 0.2)
        W = build(alg, iip, n, withZ, seed)
        calculate_step!(W, h0, nothing, nothing)
        for (kind, f) in ops
            kind === :accept ? accept_step!(W, f * W.dt, nothing, nothing) :
                reject_step!(W, f * W.dt, nothing, nothing)
        end
        return W
    end

    @testset "single rejection has the right variance" begin
        # Bridging a step of length h down to q*h must give an increment with
        # variance q*h. Getting the bridge's total-interval argument wrong instead
        # yields q*(2-q)*h.
        M = 100_000
        h = 1.0
        for alg in (:RSwM0, :RSwM1, :RSwM2, :RSwM3), q in (0.25, 0.5, 0.8)
            xs = Vector{Float64}(undef, M)
            for m in 1:M
                W = WienerProcess(
                    0.0, 0.0, nothing;
                    rswm = RSWM(adaptivealg = alg), rng = Xoshiro(m)
                )
                calculate_step!(W, h, nothing, nothing)
                reject_step!(W, q * h, nothing, nothing)
                xs[m] = W.dW
            end
            @test isapprox(var(xs) / (q * h), 1.0, atol = 6 * sqrt(2 / M))
        end
    end

    @testset "accept/reject schedules stay Brownian" begin
        M = 20_000
        for sseed in 1:2
            ops = make_schedule(seed = sseed)
            for alg in (:RSwM0, :RSwM1, :RSwM2, :RSwM3),
                    (iip, n) in ((false, 1), (true, 2)),
                    withZ in (false, true)

                ts = copy(run_path(ops, alg, iip, n, withZ, 1).t)
                nt = length(ts)
                d = iip ? n : 1
                incW = zeros(M, nt - 1, d)
                incZ = withZ ? zeros(M, nt - 1, d) : zeros(0, 0, 0)
                for m in 1:M
                    W = run_path(ops, alg, iip, n, withZ, 1000 + m)
                    @test length(W.t) == nt
                    @test maximum(abs, W.t .- ts) < 1.0e-12
                    for i in 1:(nt - 1)
                        incW[m, i, :] .= W.W[i + 1] .- W.W[i]
                        withZ && (incZ[m, i, :] .= W.Z[i + 1] .- W.Z[i])
                    end
                end
                dts = diff(ts)
                for inc in (incW, incZ)
                    isempty(inc) && continue
                    for i in 1:(nt - 1)
                        x = vec(@view inc[:, i, :])
                        @test isapprox(
                            var(x) / dts[i], 1.0,
                            atol = 6 * sqrt(2 / (M * d))
                        )
                        @test abs(mean(x)) < 6 * sqrt(dts[i] / (M * d))
                    end
                    for i in 1:(nt - 2)
                        c = cor(vec(@view inc[:, i, :]), vec(@view inc[:, i + 1, :]))
                        @test abs(c) < 6 / sqrt(M * d)
                    end
                end
            end
        end
    end

    @testset "in-place stacks never alias the scratch buffers" begin
        # ResettableStacks.push! stores by reference while copyat_or_push! copies, so
        # pushing a scratch buffer onto a stack would let later bridges silently
        # rewrite stored noise.
        ops = make_schedule(seed = 3)
        for alg in (:RSwM1, :RSwM2, :RSwM3), withZ in (false, true)
            W = build(alg, true, 3, withZ, 7)
            calculate_step!(W, 0.2, nothing, nothing)
            for (kind, f) in ops
                kind === :accept ? accept_step!(W, f * W.dt, nothing, nothing) :
                    reject_step!(W, f * W.dt, nothing, nothing)
                scratch = Any[W.dW, W.dWtilde, W.dWtmp, W.curW]
                withZ && append!(scratch, Any[W.dZ, W.dZtilde, W.dZtmp, W.curZ])
                entries = Any[]
                for S in (W.S₁, W.S₂), el in collect(S)
                    push!(entries, el[2])
                    withZ && push!(entries, el[3])
                end
                for e in entries
                    @test !any(b -> e === b, scratch)
                end
                for a in 1:length(entries), b in (a + 1):length(entries)
                    @test !(entries[a] === entries[b])
                end
            end
        end
    end

    @testset "RSwM3 S₂ decomposes the current step" begin
        # S₂ holds the breakdown of the step currently being attempted; its lengths
        # must add up to W.dt or a rejection walks back over the wrong interval.
        ops = make_schedule(seed = 4)
        for (iip, n) in ((false, 1), (true, 2))
            W = build(:RSwM3, iip, n, true, 11)
            calculate_step!(W, 0.2, nothing, nothing)
            for (kind, f) in ops
                kind === :accept ? accept_step!(W, f * W.dt, nothing, nothing) :
                    reject_step!(W, f * W.dt, nothing, nothing)
                s = sum(el[1] for el in collect(W.S₂); init = zero(W.dt))
                @test isapprox(s, W.dt, atol = 1.0e-12)
            end
        end
    end
end
