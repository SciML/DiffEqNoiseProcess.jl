function accept_step!(W::NoiseProcess,dt,setup_next=true)

  W.curW += W.dW
  W.curt += W.dt
  push!(W.W,copy(W.curW))
  push!(W.t,copy(W.curt))
  if W.Z != nothing
    W.curZ += W.dZ
    push!(W.Z,copy(W.curZ))
  end

  W.dt = dt #dtpropose
  # Setup next step
  if setup_next
    setup_next_step!(W::NoiseProcess)
  end
end

function setup_next_step!(W::NoiseProcess)
  if adaptive_alg(W)==:RSwM3
    ResettableStacks.reset!(W.S₂) #Empty W.S₂
  end
  if adaptive_alg(W)==:RSwM1
    if !isempty(W.S₁)
      W.dt,W.dW,W.dZ = pop!(W.S₁)
    else # Stack is empty
      calculate_step!(W,W.dt)
    end
  elseif adaptive_alg(W)==:RSwM2 || adaptive_alg(W)==:RSwM3
    if !(typeof(W.dW) <: AbstractArray)
      dttmp = 0.0; W.dW = 0.0
      if W.Z != nothing
        W.dZ = 0.0
      end
    else
      dttmp = 0.0; fill!(W.dW,zero(eltype(W.dW)))
      if W.Z != nothing
        fill!(W.dZ,zero(eltype(W.dZ)))
      end
    end
    while !isempty(W.S₁)
      L₁,L₂,L₃ = pop!(W.S₁)
      qtmp = (W.dt-dttmp)/L₁
      if qtmp>1
        dttmp+=L₁
        if typeof(W.dW) <: AbstractArray
          for i in eachindex(W.dW)
            W.dW[i]+=L₂[i]
            if W.Z != nothing
              W.dZ[i]+=L₃[i]
            end
          end
        else
          W.dW+=L₂
          if W.Z != nothing
            W.dZ+=L₃
          end
        end
        if adaptive_alg(W)==:RSwM3
          push!(W.S₂,(L₁,L₂,L₃))
        end
      else #Popped too far
        # Generate numbers to bridge and step perfectly
        dttmp += qtmp*L₁
        if isinplace(W)
          W.bridge(W.dWtilde,W,W.curW,L₂,qtmp,L₁)
          if W.Z != nothing
            W.bridge(W.dZtilde,W,W.curZ,L₃,qtmp,L₁)
          end
        else
          W.dWtilde = W.bridge(W,W.curW,L₂,qtmp,L₁)
          if W.Z != nothing
            W.dZtilde = W.bridge(W,W.curZ,L₃,qtmp,L₁)
          end
        end
        if typeof(W.dW) <: AbstractArray
          for i in eachindex(W.dW)
            W.dW[i] += W.dWtilde[i]
            if W.Z != nothing
              W.dZ[i] += W.dZtilde[i]
            end
          end
        else
          W.dW += W.dWtilde
          if W.Z != nothing
            W.dZ += W.dZtilde
          end
        end
        if (1-qtmp)*L₁ > W.rswm.discard_length
          if W.Z == nothing
            push!(W.S₁,((1-qtmp)*L₁,L₂-W.dWtilde,nothing))
          else
            push!(W.S₁,((1-qtmp)*L₁,L₂-W.dWtilde,L₃-W.dZtilde))
          end
          if adaptive_alg(W)==:RSwM3 && qtmp*L₁ > W.rswm.discard_length
            if W.Z == nothing
              push!(W.S₂,(qtmp*L₁,copy(W.dWtilde),nothing))
            else
              push!(W.S₂,(qtmp*L₁,copy(W.dWtilde),copy(W.dZtilde)))
            end
          end
        end
        break
      end
    end #end while empty
    dtleft = W.dt - dttmp
    if !≈(dtleft,0.0,atol=W.rswm.discard_length) #Stack emptied
      if isinplace(W)
        W.dist(W.dWtilde,W,dtleft)
        if W.Z != nothing
          W.dist(W.dZtilde,W,dtleft)
        end
      else
        W.dWtilde = W.dist(W,dtleft)
        if W.Z != nothing
          W.dZtilde = W.dist(W,dtleft)
        end
      end
      if typeof(W.dW) <: AbstractArray
        for i in eachindex(W.dW)
          W.dW[i] += W.dWtilde[i]
          if W.Z != nothing
            W.dZ[i] += W.dZtilde[i]
          end
        end
      else
        W.dW += W.dWtilde
        if W.Z != nothing
          W.dZ += W.dZtilde
        end
      end
      if adaptive_alg(W)==:RSwM3
        if W.Z == nothing
          push!(W.S₂,(dtleft,copy(W.dWtilde),nothing))
        else
          push!(W.S₂,(dtleft,copy(W.dWtilde),copy(W.dZtilde)))
        end
      end
    end
  end # End RSwM2 and RSwM3
end

function calculate_step!(W::NoiseProcess,dt)
  if isinplace(W)
    W.dist(W.dW,W,dt)
    if W.Z != nothing
      W.dist(W.dZ,W,dt)
    end
  else
    W.dW = W.dist(W,dt)
    if W.Z != nothing
      W.dZ = W.dist(W,dt)
    end
  end
  W.dt = dt
end

function reject_step!(W::NoiseProcess,dtnew)
  q = dtnew/W.dt
  if adaptive_alg(W)==:RSwM1 || adaptive_alg(W)==:RSwM2
    if isinplace(W)
      W.bridge(W.dWtilde,W,W.curW,W.curW+W.dW,q,dtnew)
      W.dWtilde .-= W.curW
      if W.Z != nothing
        W.bridge(W.dZtilde,W,W.curZ,W.curZ+W.dZ,q,dtnew)
        W.dZtilde .-= W.curZ
      end
    else
      W.dWtilde = W.bridge(W,W.curW,W.curW+W.dW,q,dtnew)-W.curW
      if W.Z != nothing
        W.dZtilde=  W.bridge(W,W.curZ,W.curZ+W.dZ,q,dtnew)-W.curZ
      end
    end
    cutLength = W.dt-dtnew
    if cutLength > W.rswm.discard_length
      if W.Z == nothing
        push!(W.S₁,(cutLength,W.dW-W.dWtilde,nothing))
      else
        push!(W.S₁,(cutLength,W.dW-W.dWtilde,W.dZ-W.dZtilde))
      end
    end
    if length(W.S₁) > W.maxstacksize
        W.maxstacksize = length(W.S₁)
    end
    if typeof(W.dW) <: AbstractArray
      copy!(W.dW,W.dWtilde)
      if W.Z!=nothing
        copy!(W.dZ,W.dZtilde)
      end
    else
      W.dW = W.dWtilde
      if W.Z != nothing
        W.dZ = W.dZtilde
      end
    end
    W.dt = dtnew
  else # RSwM3
    if !(typeof(W.dW) <: AbstractArray)
      dttmp = 0.0; W.dWtmp = 0.0
      if W.Z != nothing
        W.dZtmp = 0.0
      end
    else
      dttmp = 0.0; fill!(W.dWtmp,zero(eltype(W.dWtmp)))
      if W.Z!= nothing
        fill!(W.dZtmp,zero(eltype(W.dZtmp)))
      end
    end
    if length(W.S₂) > W.maxstacksize2
      W.maxstacksize2= length(W.S₂)
    end
    while !isempty(W.S₂)
      L₁,L₂,L₃ = pop!(W.S₂)
      if dttmp + L₁ < (1-q)*W.dt #while the backwards movement is less than chop off
        dttmp += L₁
        if typeof(W.dW) <: AbstractArray
          for i in eachindex(W.dW)
            W.dWtmp[i] += L₂[i]
            if W.Z != nothing
              W.dZtmp[i] += L₃[i]
            end
          end
        else
          W.dWtmp += L₂
          if W.Z != nothing
            W.dZtmp += L₃
          end
        end
        push!(W.S₁,(L₁,L₂,L₃))
      else
        push!(W.S₂,(L₁,L₂,L₃))
        break
      end
    end # end while
    dtK = W.dt - dttmp
    qK = q*W.dt/dtK
    if typeof(W.dW) <: AbstractArray
      for i in eachindex(W.dW)
        W.dWtmp[i] = W.dW[i] - W.dWtmp[i]
        if W.Z != nothing
          W.dZtmp[i] = W.dZ[i] - W.dZtmp[i]
        end
      end
    else
      W.dWtmp = W.dW - W.dWtmp
      if W.Z != nothing
        W.dZtmp = W.dZ - W.dZtmp
      end
    end
    if isinplace(W)
      W.bridge(W.dWtilde,W,0,W.dWtmp,qK,dtK)
      #W.dWtilde .-= W.curW
      if W.Z != nothing
        W.bridge(W.dZtilde,W,0,W.dZtmp,qK,dtK)
        #W.dZtilde .-= W.curZ
      end
    else
      W.dWtilde = W.bridge(W,0,W.dWtmp,qK,dtK)# - W.curW
      if W.Z != nothing
        W.dZtilde = W.bridge(W,0,W.dZtmp,qK,dtK)# - W.curZ
      end
    end
    cutLength = (1-qK)*dtK
    if cutLength > W.rswm.discard_length
      if W.Z == nothing
        push!(W.S₁,(cutLength,W.dWtmp-W.dWtilde,nothing))
      else
        push!(W.S₁,(cutLength,W.dWtmp-W.dWtilde,W.dZtmp-W.dZtilde))
      end
    end
    if length(W.S₁) > W.maxstacksize
        W.maxstacksize = length(W.S₁)
    end
    W.dt = dtnew
    if typeof(W.dW) <: AbstractArray
      copy!(W.dW,W.dWtilde)
      if W.Z != nothing
        copy!(W.dZ,W.dZtilde)
      end
    else
      W.dW = W.dWtilde
      if W.Z != nothing
        W.dZ = W.dZtilde
      end
    end
  end
end

function interpolate!(W::NoiseProcess,t)
  if t > W.t[end] # Steps past W
    dt = t - W.t[end]
    if isinplace(W)
      W.dist(W.dW,W,dt)
      W.curW .+= W.dW
      if W.Z != nothing
        W.dist(W.dZ,W,dt)
        W.curZ .+= W.dZ
      end
    else
      W.dW = W.dist(W,dt)
      W.curW += W.dW
      if W.Z != nothing
        W.dZ = W.dist(W,dt)
        W.curZ += W.dZ
      end
    end
    push!(W.t,t)
    out1 = copy(W.curW)
    push!(W.W,out1)
    if W.Z != nothing
      out2= copy(W.curZ)
      push!(W.Z,out2)
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
        W.bridge(new_curW,W,W0,Wh,q,h)
        if W.Z != nothing
          new_curZ = similar(W.dZ)
          W.bridge(new_curZ,W,Z0,Zh,q,h)
        else
          new_curZ = nothing
        end
      else
        new_curW = W.bridge(W,W0,Wh,q,h)
        if W.Z != nothing
          new_curZ = W.bridge(W,Z0,Zh,q,h)
        else
          new_curZ = nothing
        end
      end
      W.curW = new_curW
      insert!(W.W,i,new_curW)
      insert!(W.t,i,t)
      if W.Z != nothing
        W.curZ = new_curZ
        insert!(W.Z,i,new_curZ)
      end
      return new_curW,new_curZ
    end
  end
end

function interpolate!(out1,out2,W::NoiseProcess,t)
  if t > W.t[end] # Steps past W
    dt = t - W.t[end]
    W.dist(W.dW,W,dt)
    out1 .+= W.dW
    if W.Z != nothing
      W.dist(W.dZ,W,dt)
      out2 .+= W.dZ
    end
    push!(W.t,t)
    push!(W.W,copy(out1))
    if W.Z != nothing
      push!(W.Z,copy(out2))
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
      W.bridge(out1,W,W0,Wh,q,h)
      if W.Z != nothing
        W.bridge(out2,W,Z0,Zh,q,h)
      end
      W.curW .= out1
      insert!(W.W,i,copy(out1))
      insert!(W.t,i,t)
      if W.Z != nothing
        W.curZ .= out2
        insert!(W.Z,i,copy(out2))
      end
    end
  end
end
