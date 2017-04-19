function accept_step!(W::WienerProcess,dt)

  W.curW += W.dW
  W.curZ += W.dZ
  W.curt += W.dt
  push!(W.W,W.curW)
  push!(W.Z,W.curZ)
  push!(W.t,W.curt)

  W.dt = dt #dtpropose

  #modify_dt_for_tstops!(W)
  # End with #W.sqdt = sqrt(abs(W.dt))


  # Setup next step
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
      dttmp = 0.0; W.dW = 0.0; W.dZ = 0.0
    else
      dttmp = 0.0; fill!(W.dW,zero(eltype(W.dW))); fill!(W.dZ,zero(eltype(W.dZ)))
    end
    while !isempty(W.S₁)
      L₁,L₂,L₃ = pop!(W.S₁)
      qtmp = (W.dt-dttmp)/L₁
      if qtmp>1
        dttmp+=L₁
        if typeof(W.dW) <: AbstractArray
          for i in eachindex(W.dW)
            W.dW[i]+=L₂[i]; W.dZ[i]+=L₃[i]
          end
        else
          W.dW+=L₂; W.dZ+=L₃
        end
        if adaptive_alg(W)==:RSwM3
          push!(W.S₂,(L₁,L₂,L₃))
        end
      else #Popped too far
        W.dWtilde = W.bridge(W.curW,L₂,qtmp,L₁)-W.curW
        W.dZtilde = W.bridge(W.curZ,L₃,qtmp,L₁)-W.curZ
        if typeof(W.dW) <: AbstractArray
          for i in eachindex(W.dW)
            W.dW[i] += W.dWtilde[i]; W.ΔZ[i] += W.dZtilde[i]
          end
        else
          W.dW += W.dWtilde; W.dZ += W.dZtilde
        end
        if (1-qtmp)*L₁ > W.rswm.discard_length
          push!(W.S₁,((1-qtmp)*L₁,L₂-W.dWtilde,L₃-W.dZtilde))
          if adaptive_alg(W)==:RSwM3 && qtmp*L₁ > W.rswm.discard_length
            push!(W.S₂,(qtmp*L₁,copy(W.dWtilde),copy(W.dZtilde)))
          end
        end
        break
      end
    end #end while empty
    dtleft = W.dt - dttmp
    if dtleft != 0 #Stack emptied
      W.dWtilde = W.dist(W,dtleft)
      W.dZtilde = W.dist(W,dtleft)
      if typeof(W.dW) <: AbstractArray
        for i in eachindex(W.dW)
          W.dW[i] += W.dWtilde[i]; W.dZ[i] += W.dZtilde[i]
        end
      else
        W.dW += W.dWtilde; W.dZ += W.dZtilde
      end
      if adaptive_alg(W)==:RSwM3
        push!(W.S₂,(dtleft,copy(W.dWtilde),copy(W.dZtilde)))
      end
    end
  end # End RSwM2 and RSwM3
end

function calculate_step!(W::WienerProcess,dt)
  W.dW = W.dist(W,dt)
  W.dZ = W.dist(W,dt)
  W.dt = dt
end

function reject_step!(W::WienerProcess,dtnew)
  q = dtnew/W.dt
  if adaptive_alg(W)==:RSwM1 || adaptive_alg(W)==:RSwM2
    W.dWtilde = W.bridge(W.curW,W.curW+W.dW,q,dtnew)-W.curW
    W.dZtilde=  W.bridge(W.curZ,W.curZ+W.dZ,q,dtnew)-W.curZ
    cutLength = W.dt-dtnew
    if cutLength > W.rswm.discard_length
      push!(W.S₁,(cutLength,W.dW-W.dWtilde,W.dZ-W.dZtilde))
    end
    if length(W.S₁) > W.maxstacksize
        W.maxstacksize = length(W.S₁)
    end
    if typeof(W.dW) <: AbstractArray
      copy!(W.dW,W.dWtilde); copy!(W.dZ,W.dZtilde)
    else
      W.dW = W.dWtilde; W.dZ = W.dZtilde
    end
    W.dt = dtnew
  else # RSwM3
    if !(typeof(W.dW) <: AbstractArray)
      dttmp = 0.0; W.dWtmp = 0.0; W.dZtmp = 0.0
    else
      dttmp = 0.0; fill!(W.dWtmp,zero(eltype(W.dWtmp))); fill!(W.dZtmp,zero(eltype(W.dZtmp)))
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
            W.dWtmp[i] += L₂[i]; W.dZtmp[i] += L₃[i]
          end
        else
          W.dWtmp += L₂; W.dZtmp += L₃
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
        W.dWtmp[i] = W.dW[i] - W.dWtmp[i]; W.dZtmp[i] = W.dZ[i] - W.dZtmp[i]
      end
    else
      W.dWtmp = W.dW - W.dWtmp; W.dZtmp = W.dZ - W.dZtmp
    end
    W.dWtilde = W.bridge(0,W.dWtmp,qK,dtK)
    W.dZtilde = W.bridge(0,W.dZtmp,qK,dtK)
    cutLength = (1-qK)*dtK
    if cutLength > W.rswm.discard_length
      push!(W.S₁,(cutLength,W.dWtmp-W.dWtilde,W.dZtmp-W.dZtilde))
    end
    if length(W.S₁) > W.maxstacksize
        W.maxstacksize = length(W.S₁)
    end
    W.dt = dtnew
    if typeof(W.dW) <: AbstractArray
      copy!(W.dW,W.dWtilde); copy!(W.dZ,W.dZtilde)
    else
      W.dW = W.dWtilde;  W.dZ = W.dZtilde
    end
  end

  # W.sqdt = sqrt(abs(W.dt))
end

function interpolate!(W::WienerProcess,t)
  if t > W.t[end] # Steps past W
    error("Cannot extrapolate with the interpolation")
  else # Bridge
    i = searchsortedfirst(W.t,t)
    if t == W.t[i]
      return W.W[i]
    else
      W0,Wh = W.W[i-1],W.W[i]
      Z0,Zh = W.Z[i-1],W.Z[i]
      h = W.t[i]-W.t[i-1]
      q = (t-W.t[i-1])/h
      new_curW = W.bridge(W0,Wh,q,h)
      new_curZ = W.bridge(Z0,Zh,q,h)
      W.curW = new_curW
      W.curZ = new_curZ
      insert!(W.W,i,new_curW)
      insert!(W.Z,i,new_curZ)
      insert!(W.t,i,t)
      return new_curW
    end
  end
end
