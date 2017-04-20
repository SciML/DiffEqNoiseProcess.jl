function interpolate!(W::NoiseWrapper,t)
  W.source(t)
end

function interpolate!(out1,out2,W::NoiseWrapper,t)
  W.source(out1,out2,t)
end

function calculate_step!(W::NoiseWrapper,dt)
  if isinplace(W)
    W(W.dW,W.dZ,W.t+dt)
    W.dW .-= W.curW
    W.dZ .-= W.curZ
  else
    new_W, new_Z = W(W.curt+dt)
    W.dW = new_W - W.curW
    if W.Z != nothing
      W.dZ = new_Z - W.curZ
    end
  end
  W.dt = dt
end

Base.getindex(W::NoiseWrapper,i::Int) = W.W[i]

function accept_step!(W::NoiseWrapper,dt,setup_next=true)

  W.curW += W.dW
  W.curt += W.dt
  push!(W.W,W.curW)
  push!(W.t,W.curt)
  if W.Z != nothing
    W.curZ += W.dZ
    push!(W.Z,W.curZ)
  end

  W.dt = dt #dtpropose
  if setup_next
    calculate_step!(W,dt)
  end
end

function reject_step!(W::NoiseWrapper,dtnew)
  calculate_step!(W::NoiseWrapper,dtnew)
end
