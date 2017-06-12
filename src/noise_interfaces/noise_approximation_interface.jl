function save_noise!(W::NoiseApproximation)

end

function interpolate!(W::NoiseApproximation,t)
  W.Z != nothing ? z = W.source2.sol(t) : z = nothing
  W.source1.sol(t),z
end

function interpolate!(out1,out2,W::NoiseApproximation,t)
  W.source1.sol(out1,t), W.source2.sol(out2,t)
end

function calculate_step!(W::NoiseApproximation,dt)
  W.dt = dt
  t = W.curt+dt

  add_tstop!(W.source1,t)
  step!(W.source1)
  if W.Z != nothing
    add_tstop!(W.source2,t)
    step!(W.source2)
  end

  if isinplace(W)
    W.dW .= W.source1.u .- W.curW
    if W.Z != nothing
      W.dZ .= W.source2.u .- W.curZ
    end
  else
    W.dW = W.source1.u - W.curW
    if W.Z != nothing
      W.dZ = W.source2.u - W.curZ
    end
  end
end

function accept_step!(W::NoiseApproximation,dt,setup_next=true)
  if isinplace(W)
    W.curW .+= W.dW
  else
    W.curW += W.dW
  end
  W.curt += W.dt
  if W.Z != nothing
    if isinplace(W)
      W.curZ .+= W.dZ
    else
      W.curZ += W.dZ
    end
  end

  W.dt = dt #dtpropose
  if setup_next
    calculate_step!(W,dt)
  end
end

function reject_step!(W::NoiseApproximation,dtnew)
  W.dt = dtnew
  if isinplace(W)
    W.source1(W.dW,W.curt+dtnew)
    W.dW .-= W.curW
    if W.Z != nothing
      W.source2(W.dZ,W.curt+dtnew)
      W.dZ .-= W.curZ
    end
  else
    W.dW = W.source1(W.curt+dtnew) - W.curW
    if W.Z != nothing
      W.dZ = W.source2(W.curt+dtnew) - W.curZ
    end
  end
end

function setup_next_step!(W::NoiseApproximation)
  calculate_step!(W,W.dt)
end
