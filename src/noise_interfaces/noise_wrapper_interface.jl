function save_noise!(W::NoiseWrapper)

end

function interpolate!(W::NoiseWrapper,t)
  W.source(t)
end

function interpolate!(out1,out2,W::NoiseWrapper,t)
  W.source(out1,out2,t)
end

function calculate_step!(W::NoiseWrapper,dt)
  if isinplace(W)
    W(W.dW,W.dZ,W.curt+dt)
    W.dW .-= W.curW
    if W.Z != nothing
      W.dZ .-= W.curZ
    end
  else
    new_W, new_Z = W(W.curt+dt)
    W.dW = new_W - W.curW
    if W.Z != nothing
      W.dZ = new_Z - W.curZ
    end
  end
  W.dt = dt
end

function accept_step!(W::NoiseWrapper,dt,setup_next=true)
  if isinplace(W)
    W.curW .+= W.dW
  else
    W.curW += W.dW
  end
  W.curt += W.dt
  push!(W.W,copy(W.curW))
  push!(W.t,copy(W.curt))
  if W.Z != nothing
    if isinplace(W)
      W.curZ .+= W.dZ
    else
      W.curZ += W.dZ
    end

    push!(W.Z,copy(W.curZ))
  end

  W.dt = dt #dtpropose
  if setup_next
    calculate_step!(W,dt)
  end
end

function reject_step!(W::NoiseWrapper,dtnew)
  calculate_step!(W::NoiseWrapper,dtnew)
end

function setup_next_step!(W::NoiseWrapper)
  calculate_step!(W,W.dt)
end
