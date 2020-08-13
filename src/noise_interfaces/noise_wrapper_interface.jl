function save_noise!(W::NoiseWrapper)

end

function interpolate!(W::NoiseWrapper,u,p,t; reverse=false)
  if reverse
    interpolate!(W,u,p,t,reverse=reverse)
  else
    W.source(u,p,t)
  end
end

function interpolate!(out1,out2,W::NoiseWrapper,u,p,t; reverse=false)
  if reverse
    interpolate!(out1,out2,W.source,u,p,t,reverse=reverse)
  else
    W.source(out1,out2,u,p,t)
  end
end

function calculate_step!(W::NoiseWrapper,dt,u,p)
  if isinplace(W)
    W(W.dW,W.dZ,u,p,W.curt+dt)
    W.dW .-= W.curW
    if W.Z != nothing
      W.dZ .-= W.curZ
    end
  else
    new_W, new_Z = W(u,p,W.curt+dt)
    W.dW = new_W - W.curW
    if W.Z != nothing
      W.dZ = new_Z - W.curZ
    end
  end
  W.dt = dt
  return nothing
end

function accept_step!(W::NoiseWrapper,dt,u,p,setup_next=true)
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
    calculate_step!(W,dt,u,p)
  end
  return nothing
end

function reject_step!(W::NoiseWrapper,dtnew,u,p)
  calculate_step!(W::NoiseWrapper,dtnew,u,p)
  return nothing
end

function setup_next_step!(W::NoiseWrapper,u,p)
  calculate_step!(W,W.dt,u,p)
  return nothing
end
