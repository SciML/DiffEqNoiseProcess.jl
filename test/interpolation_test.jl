@testset "Noise Interpolation Test" begin

    using DiffEqNoiseProcess

    W = WienerProcess(0.0, 0.0, 0.0)

    dt = 0.1
    calculate_step!(W, dt, nothing, nothing)

    for i = 1:10
        accept_step!(W, dt, nothing, nothing)
    end

    dt = dt / 100
    for t = dt:dt:1-dt
        W(t)
    end

    dt = dt / 100
    for t = dt:dt:1-dt
        W(t)
    end

end
