@inline function save_noise!(W::SimpleNoiseProcess)
  if W.t != W.curt
    push!(W.W,copy(W.curW))
    push!(W.t,copy(W.curt))
    if W.Z != nothing
      push!(W.Z,copy(W.curZ))
    end
  end
end

@inline function accept_step!(W::SimpleNoiseProcess,dt,u,p,setup_next=true)

  W.curt += W.dt
  W.iter += 1

  if isinplace(W)
    @.. W.curW += W.dW
    if W.Z != nothing
      @.. W.curZ += W.dZ
    end
  else
    W.curW += W.dW
    if W.Z != nothing
      W.curZ += W.dZ
    end
  end

  if W.save_everystep
    push!(W.W,copy(W.curW))
    push!(W.t,copy(W.curt))
    if W.Z != nothing
      push!(W.Z,copy(W.curZ))
    end
  end

  W.dt = dt #dtpropose
  # Setup next step
  if setup_next
    setup_next_step!(W::SimpleNoiseProcess,u,p)
  end
end

@inline function setup_next_step!(W::SimpleNoiseProcess,u,p)
  calculate_step!(W,W.dt,u,p)
end

@inline function calculate_step!(W::SimpleNoiseProcess,dt,u,p)
  if isinplace(W)
    W.dist(W.dW,W,dt,u,p,W.curt,W.rng)
    if W.Z != nothing
      W.dist(W.dZ,W,dt,u,p,W.curt,W.rng)
    end
  else
    W.dW = W.dist(W,dt,u,p,W.curt,W.rng)
    if W.Z != nothing
      W.dZ = W.dist(W,dt,u,p,W.curt,W.rng)
    end
  end
  W.dt = dt
end

@inline function reject_step!(W::SimpleNoiseProcess,dtnew,u,p)
  error("SimpleNoiseProcess cannot be used with adaptivity rejections")
end

@inline function interpolate!(W::SimpleNoiseProcess,t)
  if sign(W.dt)*t > sign(W.dt)*W.t[end] # Steps past W
    dt = t - W.t[end]
    if isinplace(W)
      W.dist(W.dW,W,dt,u,p,t,W.rng)
      W.curW .+= W.dW
      if W.Z != nothing
        W.dist(W.dZ,W,dt,u,p,t,W.rng)
        W.curZ .+= W.dZ
      end
    else
      W.dW = W.dist(W,dt,u,p,t,W.rng)
      W.curW += W.dW
      if W.Z != nothing
        W.dZ = W.dist(W,dt,u,p,t,W.rng)
        W.curZ += W.dZ
      end
    end
    out1 = copy(W.curW)
    if W.save_everystep
      push!(W.t,t)
      push!(W.W,out1)
    end
    if W.Z != nothing
      out2= copy(W.curZ)
      if W.save_everystep
        push!(W.Z,out2)
      end
    else
      out2 = nothing
    end
    return out1,out2
  else # Bridge
    error("SimpleNoiseProcess cannot interpolate")
  end
end

@inline function interpolate!(out1,out2,W::SimpleNoiseProcess,t)
  if sign(W.dt)*t > sign(W.dt)*W.t[end] # Steps past W
    dt = t - W.t[end]
    W.dist(W.dW,W,dt,u,p,t,W.rng)
    out1 .+= W.dW
    if W.Z != nothing
      W.dist(W.dZ,W,dt,u,p,t,W.rng)
      out2 .+= W.dZ
    end
    if W.save_everystep
      push!(W.t,t)
      push!(W.W,copy(out1))
      if W.Z != nothing
        push!(W.Z,copy(out2))
      end
    end
  else # Bridge
    error("SimpleNoiseProcess cannot interpolate")
  end
end
