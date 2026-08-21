@setup_workload begin
    process = WienerProcess(0.0, 0.0; rng = Xoshiro(0))

    @compile_workload begin
        calculate_step!(process, 0.1, nothing, nothing)
        accept_step!(process, 0.1, nothing, nothing, false)
        process(0.05)
        reject_step!(process, 0.05, nothing, nothing)
        accept_step!(process, 0.05, nothing, nothing, false)
        process(0.1)
    end
end
