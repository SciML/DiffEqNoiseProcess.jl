@inline function save_noise!(W::SimpleNoiseProcess)
  if W.t != W.curt
    push!(W.W,copy(W.curW))
    push!(W.t,copy(W.curt))
    if W.Z != nothing
      push!(W.Z,copy(W.curZ))
    end
  end
  return nothing
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
  return nothing
end

@inline function setup_next_step!(W::SimpleNoiseProcess,u,p)
  calculate_step!(W,W.dt,u,p)
  return nothing
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
  return nothing
end

@inline function reject_step!(W::SimpleNoiseProcess,dtnew,u,p)
  error("SimpleNoiseProcess cannot be used with adaptivity rejections")
end

@inline function interpolate!(W::SimpleNoiseProcess,u,p,t)
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
    i = searchsortedfirst(W.t,t)
    if t == W.t[i]
      if isinplace(W)
        W.curW .= W.W[i]
      else
        W.curW = W.W[i]
      end
      if W.Z != nothing
        if isinplace(W)
          W.curZ .= W.Z[i]
        else
          W.curZ = W.Z[i]
        end
        return copy(W.curW),copy(W.curZ)
      else
        return copy(W.curW),nothing
      end
    else
      W0,Wh = W.W[i-1],W.W[i]
      if W.Z != nothing
        Z0,Zh = W.Z[i-1],W.Z[i]
      end
      h = W.t[i]-W.t[i-1]
      q = (t-W.t[i-1])/h
      if isinplace(W)
        new_curW = similar(W.dW)
        W.bridge(new_curW,W,W0,Wh,q,h,u,p,t,W.rng)
        if iscontinuous(W)
          @. new_curW += (1-q)*W0
        else
          @. new_curW += W0
        end
        if W.Z != nothing
          new_curZ = similar(W.dZ)
          W.bridge(new_curZ,W,Z0,Zh,q,h,u,p,t,W.rng)
          if iscontinuous(W)
            @. new_curZ += (1-q)*Z0
          else
            @. new_curZ += Z0
          end
        else
          new_curZ = nothing
        end
      else
        new_curW = W.bridge(W,W0,Wh,q,h,u,p,t,W.rng)
        if iscontinuous(W)
          # This should actually be based on the function for computing the mean
          # flow of the noise process, but for now we'll just handle Wiener and
          # Poisson
          new_curW += (1-q)*W0
        else
          new_curW += W0
        end
        if W.Z != nothing
          new_curZ = W.bridge(W,Z0,Zh,q,h,u,p,t,W.rng)
          if iscontinuous(W)
            new_curZ += (1-q)*Z0
          else
            new_curZ += Z0
          end
        else
          new_curZ = nothing
        end
      end
      W.curW = new_curW
      if W.save_everystep
        insert!(W.W,i,new_curW)
        insert!(W.t,i,t)
      end
      if W.Z != nothing
        W.curZ = new_curZ
        if W.save_everystep
          insert!(W.Z,i,new_curZ)
        end
      end
      return new_curW,new_curZ
    end
  end
end

@inline function interpolate!(out1,out2,W::SimpleNoiseProcess,u,p,t)
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
    i = searchsortedfirst(W.t,t)
    if t == W.t[i]
      out1 .= W.W[i]
      if W.Z != nothing
        out2 .= W.Z[i]
      end
    else
      W0,Wh = W.W[i-1],W.W[i]
      if W.Z != nothing
        Z0,Zh = W.Z[i-1],W.Z[i]
      end
      h = W.t[i]-W.t[i-1]
      q = (t-W.t[i-1])/h
      W.bridge(out1,W,W0,Wh,q,h,u,p,t,W.rng)
      out1 .+= (1-q)*W0
      if W.Z != nothing
        W.bridge(out2,W,Z0,Zh,q,h,u,p,t,W.rng)
        out2 .+= (1-q)*Z0
      end
      W.curW .= out1
      if W.save_everystep
        insert!(W.W,i,copy(out1))
        insert!(W.t,i,t)
      end
      if W.Z != nothing
        W.curZ .= out2
        if W.save_everystep
          insert!(W.Z,i,copy(out2))
        end
      end
    end
  end
end
