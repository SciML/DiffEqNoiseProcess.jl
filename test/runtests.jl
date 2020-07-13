using Test

@time begin
  include("interpolation_test.jl")
  include("RSwM1_test.jl")
  include("RSwM2_test.jl")
  include("RSwM3_test.jl")
  include("noise_wrapper.jl")
  include("noise_function.jl")
  include("noise_grid.jl")
  include("noise_approximation.jl")
  include("sde_noise_wrapper.jl")
  include("multi_dim.jl")
  include("geometric_bm.jl")
  include("compoundpoisson.jl")
  include("ornstein.jl")
  include("bridge_test.jl")
  include("sde_adaptivedistribution_tests.jl")
  include("reversal_test.jl")
  include("two_processes.jl")
  include("extraction_test.jl")
  include("restart_test.jl")
end
