function save_noise!(W::NoiseFunction)

end

Base.show(io::IO, A::NoiseFunction) =
           invoke(show, Tuple{typeof(io), Any}, io, A)
Base.show(io::IO, ::MIME"text/plain", A::NoiseFunction) = show(io, A)

function interpolate!(W::NoiseFunction,u,p,t)
  W.W(u,p,t)
end

function interpolate!(out1,out2,W::NoiseFunction,u,p,t)
  W.W(out1,out2,u,p,t)
end

function calculate_step!(W::NoiseFunction,dt,u,p)
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

function accept_step!(W::NoiseFunction,dt,u,p,setup_next=true)
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
    calculate_step!(W,dt,u,p)
  end
  return nothing
end

function reject_step!(W::NoiseFunction,dtnew,u,p)
  calculate_step!(W,dtnew,u,p)
  return nothing
end

function setup_next_step!(W::NoiseFunction,u,p)
  calculate_step!(W,W.dt,u,p)
  return nothing
end
