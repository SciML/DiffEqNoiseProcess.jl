solve(prob::AbstractNoiseProblem,::Void,args...;dt=0.0,kwargs...) = solve(prob,args...;dt=dt,kwargs...)

function solve(prob::AbstractNoiseProblem,args...;dt=0.0,kwargs...)
  if dt == 0.0 || dt == nothing
    error("dt must be provided to simulate a noise process. Please pass dt=...")
  end
  W = deepcopy(prob.noise)
  W.curt = prob.tspan[1]
  calculate_step!(W,dt)
  tType = typeof(W.curt)
  while W.curt < prob.tspan[2]
    if tType <: AbstractFloat && abs(W.curt + dt - prob.tspan[2]) < 100*eps(dt) # Correct the end due to floating point error
      dt = prob.tspan[2]- W.curt
      accept_step!(W,dt)
      W.curt = prob.tspan[2]
    else
      accept_step!(W,dt)
    end
  end
  W
end
