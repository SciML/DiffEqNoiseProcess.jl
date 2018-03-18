DiffEqBase.has_reinit(i::AbstractNoiseProcess) = true
function DiffEqBase.reinit!(W::Union{NoiseProcess,NoiseApproximation},dt;
                            t0 = W.t[1],
                            erase_sol = true,
                            setup_next = false)

  if erase_sol
    resize!(W.t,1)
    resize!(W.W,1)
    if W.Z != nothing
      resize!(W.Z,1)
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
end

function DiffEqBase.reinit!(W::AbstractNoiseProcess,dt;
                            t0 = W.t[1],
                            erase_sol = true,
                            setup_next = false)
  W.curt = t0
  W.dt = dt
  if typeof(W) <: NoiseGrid
    if isinplace(W)
      W.curW .= W.W[1]
      if W.Z != nothing
        W.curZ .= W.Z[1]
      end
    else
      W.curW = W.W[1]
      if W.Z != nothing
        W.curZ = W.Z[1]
      end
    end
    W.step_setup = true
  end
  setup_next && setup_next_step!(W)
end
