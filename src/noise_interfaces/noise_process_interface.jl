@inline function save_noise!(W::NoiseProcess)
  if W.t != W.curt
    push!(W.W,copy(W.curW))
    push!(W.t,copy(W.curt))
    if W.Z != nothing
      push!(W.Z,copy(W.curZ))
    end
  end
  return nothing
end

@inline function accept_step!(W::NoiseProcess,dt,u,p,setup_next=true)

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
    setup_next_step!(W::NoiseProcess,u,p)
  end
  return nothing
end

@inline function setup_next_step!(W::NoiseProcess,u,p)
  if adaptive_alg(W)==:RSwM3
    ResettableStacks.reset!(W.S₂) #Empty W.S₂
  end
  if adaptive_alg(W)==:RSwM1
    if !isempty(W.S₁)
      W.dt,W.dW,W.dZ = pop!(W.S₁)
    else # Stack is empty
      calculate_step!(W,W.dt,u,p)
    end
  elseif adaptive_alg(W)==:RSwM2 || adaptive_alg(W)==:RSwM3
    if !isinplace(W)
      dttmp = zero(W.dt); W.dW = zero(W.dW)
      if W.Z != nothing
        W.dZ = zero(W.dZ)
      end
    else
      dttmp = zero(W.dt); fill!(W.dW,zero(eltype(W.dW)))
      if W.Z != nothing
        fill!(W.dZ,zero(eltype(W.dZ)))
      end
    end
    while !isempty(W.S₁)
      L₁,L₂,L₃ = pop!(W.S₁)
      qtmp = (W.dt-dttmp)/L₁
      if qtmp>1
        dttmp+=L₁
        if isinplace(W)
          @.. W.dW+=L₂
          if W.Z != nothing
            @.. W.dZ+=L₃
          end
        else
          W.dW+=L₂
          if W.Z != nothing
            W.dZ+=L₃
          end
        end
        if adaptive_alg(W)==:RSwM3
          if L₁ > W.rswm.discard_length
            push!(W.S₂,(L₁,L₂,L₃))
          end
        end
      else #Popped too far
        # Generate numbers to bridge and step perfectly
        dttmp += qtmp*L₁
        if isinplace(W)
          W.bridge(W.dWtilde,W,W.curW,L₂,qtmp,L₁,u,p,W.curt,W.rng)
          if W.Z != nothing
            W.bridge(W.dZtilde,W,W.curZ,L₃,qtmp,L₁,u,p,W.curt,W.rng)
          end
        else
          W.dWtilde = W.bridge(W.dW,W,W.curW,L₂,qtmp,L₁,u,p,W.curt,W.rng)
          if W.Z != nothing
            W.dZtilde = W.bridge(W.dZ,W,W.curZ,L₃,qtmp,L₁,u,p,W.curt,W.rng)
          end
        end
        if isinplace(W)
          @.. W.dW += W.dWtilde
          if W.Z != nothing
            @.. W.dZ += W.dZtilde
          end
        else
          W.dW += W.dWtilde
          if W.Z != nothing
            W.dZ += W.dZtilde
          end
        end
        if (1-qtmp)*L₁ > W.rswm.discard_length
          if isinplace(W)
            @.. L₂ -= W.dWtilde
            if W.Z != nothing
              @.. L₃ -= W.dZtilde
            end
          else
            L₂ -= W.dWtilde
            if W.Z != nothing
              L₃ -= W.dZtilde
            end
          end
          if W.Z == nothing
            push!(W.S₁,((1-qtmp)*L₁,L₂,nothing))
          else
            push!(W.S₁,((1-qtmp)*L₁,L₂,L₃))
          end
          if adaptive_alg(W)==:RSwM3 && qtmp*L₁ > W.rswm.discard_length
            if W.Z == nothing
              ResettableStacks.copyat_or_push!(W.S₂,(qtmp*L₁,W.dWtilde,nothing))
            else
              ResettableStacks.copyat_or_push!(W.S₂,(qtmp*L₁,W.dWtilde,W.dZtilde))
            end
          end
        end
        break
      end
    end #end while empty
    # This is a control variable so do not diff through it
    dtleft = DiffEqBase.ODE_DEFAULT_NORM(W.dt - dttmp,W.curt)
    if dtleft > W.rswm.discard_length #Stack emptied
      if isinplace(W)
        W.dist(W.dWtilde,W,dtleft,u,p,W.curt,W.rng)
        if W.Z != nothing
          W.dist(W.dZtilde,W,dtleft,u,p,W.curt,W.rng)
        end
      else
        W.dWtilde = W.dist(W.dW,W,dtleft,u,p,W.curt,W.rng)
        if W.Z != nothing
          W.dZtilde = W.dist(W.dZ,W,dtleft,u,p,W.curt,W.rng)
        end
      end
      if isinplace(W)
        @.. W.dW += W.dWtilde
        if W.Z != nothing
          @.. W.dZ += W.dZtilde
        end
      else
        W.dW += W.dWtilde
        if W.Z != nothing
          W.dZ += W.dZtilde
        end
      end
      if adaptive_alg(W)==:RSwM3
        if W.Z == nothing
          ResettableStacks.copyat_or_push!(W.S₂,(dtleft,W.dWtilde,nothing))
        else
          ResettableStacks.copyat_or_push!(W.S₂,(dtleft,W.dWtilde,W.dZtilde))
        end
      end
    end
  end # End RSwM2 and RSwM3
  return nothing
end

@inline function calculate_step!(W::NoiseProcess,dt,u,p)
  if isinplace(W)
    W.dist(W.dW,W,dt,u,p,W.curt,W.rng)
    if W.Z != nothing
      W.dist(W.dZ,W,dt,u,p,W.curt,W.rng)
    end
  else
    W.dW = W.dist(W.dW,W,dt,u,p,W.curt,W.rng)
    if W.Z != nothing
      W.dZ = W.dist(W.dZ,W,dt,u,p,W.curt,W.rng)
    end
  end
  W.dt = dt
  return nothing
end

@inline function reject_step!(W::NoiseProcess,dtnew,u,p)
  q = dtnew/W.dt
  if adaptive_alg(W)==:RSwM1 || adaptive_alg(W)==:RSwM2
    if isinplace(W)
      W.bridge(W.dWtilde,W,0,W.dW,q,dtnew,u,p,W.curt,W.rng)
      if W.Z != nothing
        W.bridge(W.dZtilde,W,0,W.dZ,q,dtnew,u,p,W.curt,W.rng)
      end
    else
      W.dWtilde = W.bridge(W.dW,W,0,W.dW,q,dtnew,u,p,W.curt,W.rng)
      if W.Z != nothing
        W.dZtilde=  W.bridge(W.dZ,W,0,W.dZ,q,dtnew,u,p,W.curt,W.rng)
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
    if isinplace(W)
      copyto!(W.dW,W.dWtilde)
      if W.Z!=nothing
        copyto!(W.dZ,W.dZtilde)
      end
    else
      W.dW = W.dWtilde
      if W.Z != nothing
        W.dZ = W.dZtilde
      end
    end
    W.dt = dtnew
  else # RSwM3
    if !isinplace(W)
      dttmp = zero(W.dt); W.dWtmp = zero(W.dW)
      if W.Z != nothing
        W.dZtmp = zero(W.dZtmp)
      end
    else
      dttmp = zero(W.dt); fill!(W.dWtmp,zero(eltype(W.dWtmp)))
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
        if isinplace(W)
          @.. W.dWtmp += L₂
          if W.Z != nothing
            @.. W.dZtmp += L₃
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
    if isinplace(W)
      @.. W.dWtmp = W.dW - W.dWtmp
      if W.Z != nothing
        @.. W.dZtmp = W.dZ - W.dZtmp
      end
    else
      W.dWtmp = W.dW - W.dWtmp
      if W.Z != nothing
        W.dZtmp = W.dZ - W.dZtmp
      end
    end
    if isinplace(W)
      W.bridge(W.dWtilde,W,0,W.dWtmp,qK,dtK,u,p,W.curt,W.rng)
      #W.dWtilde .-= W.curW
      if W.Z != nothing
        W.bridge(W.dZtilde,W,0,W.dZtmp,qK,dtK,u,p,W.curt,W.rng)
        #W.dZtilde .-= W.curZ
      end
    else
      W.dWtilde = W.bridge(W.dW,W,0,W.dWtmp,qK,dtK,u,p,W.curt,W.rng)# - W.curW
      if W.Z != nothing
        W.dZtilde = W.bridge(W.dZ,W,0,W.dZtmp,qK,dtK,u,p,W.curt,W.rng)# - W.curZ
      end
    end
    # This is a control variable so do not diff through it
    cutLength = DiffEqBase.ODE_DEFAULT_NORM((1-qK)*dtK,W.curt)
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
    if isinplace(W)
      copyto!(W.dW,W.dWtilde)
      if W.Z != nothing
        copyto!(W.dZ,W.dZtilde)
      end
    else
      W.dW = W.dWtilde
      if W.Z != nothing
        W.dZ = W.dZtilde
      end
    end
  end
  return nothing
end

@inline function interpolate!(W::NoiseProcess,u,p,t; reverse=false)
  if sign(W.dt)*t > sign(W.dt)*W.t[end] || (sign(W.dt)*t < sign(W.dt)*W.t[1] && reverse) # Steps past W (forward time || backward time)
    if reverse
      dt = t - W.t[1]
    else
      dt = t - W.t[end]
    end
    if isinplace(W)
      W.dist(W.dW,W,dt,u,p,t,W.rng)
      W.curW .+= W.dW
      if W.Z != nothing
        W.dist(W.dZ,W,dt,u,p,t,W.rng)
        W.curZ .+= W.dZ
      end
    else
      W.dW = W.dist(W.dW,W,dt,u,p,t,W.rng)
      W.curW += W.dW
      if W.Z != nothing
        W.dZ = W.dist(W.dZ,W,dt,u,p,t,W.rng)
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
    if t isa Union{Rational,Integer}
      i = searchsortedfirst(W.t,t)
    else
      i = searchsortedfirst(W.t,t-eps(typeof(t)))
    end
    if t isa Union{Rational,Integer} && t ==  W.t[i]
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
    elseif !(t isa Union{Rational,Integer}) && isapprox(t, W.t[i]; atol = 100eps(typeof(t)), rtol = 100eps(t))
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
      if reverse
        W0,Wh = W.W[i],W.W[i-1]
        if W.Z != nothing
          Z0,Zh = W.Z[i+1],W.Z[i]
        end
        h = W.t[i-1]-W.t[i]
        q = (t-W.t[i])/h
      else
        W0,Wh = W.W[i-1],W.W[i]
        if W.Z != nothing
          Z0,Zh = W.Z[i-1],W.Z[i]
        end
        h = W.t[i]-W.t[i-1]
        q = (t-W.t[i-1])/h
      end

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
        new_curW = W.bridge(W.dW,W,W0,Wh,q,h,u,p,t,W.rng)
        if iscontinuous(W)
          # This should actually be based on the function for computing the mean
          # flow of the noise process, but for now we'll just handle Wiener and
          # Poisson
          new_curW += (1-q)*W0
        else
          new_curW += W0
        end
        if W.Z != nothing
          new_curZ = W.bridge(W.dZ,W,Z0,Zh,q,h,u,p,t,W.rng)
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

@inline function interpolate!(out1,out2,W::NoiseProcess,u,p,t; reverse=false)
  if (sign(W.dt)*t > sign(W.dt)*W.t[end] && !reverse) || (sign(W.dt)*t < sign(W.dt)*W.t[1] && reverse) # Steps past W (forward time || backward time)
    if reverse
      dt = t - W.t[1]
    else
      dt = t - W.t[end]
    end
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
    # if abs(dt) < 1e-14
    #   @. out1 = W.curW + W.dW
    #   if W.Z != nothing
    #     @. out2 = W.curZ + W.dZ
    #   end
    # end

  else # Bridge
    if t isa Union{Rational,Integer}
      i = searchsortedfirst(W.t,t)
    else
      i = searchsortedfirst(W.t,t-eps(typeof(t)))
    end

    if t isa Union{Rational,Integer} && t ==  W.t[i]
      out1 .= W.W[i]
      if W.Z != nothing
        out2 .= W.Z[i]
      end
    elseif !(t isa Union{Rational,Integer}) && isapprox(t, W.t[i]; atol = 100eps(typeof(t)), rtol = 100eps(t))
      out1 .= W.W[i]
      if W.Z != nothing
        out2 .= W.Z[i]
      end
    else
      if reverse
        W0,Wh = W.W[i],W.W[i-1]
        if W.Z != nothing
          Z0,Zh = W.Z[i+1],W.Z[i]
        end
        h = W.t[i-1]-W.t[i]
        q = (t-W.t[i])/h
      else
        W0,Wh = W.W[i-1],W.W[i]
        if W.Z != nothing
          Z0,Zh = W.Z[i-1],W.Z[i]
        end
        h = W.t[i]-W.t[i-1]
        q = (t-W.t[i-1])/h
      end

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
  return nothing
end

function resize_stack!(W::NoiseProcess,i)
  for j in eachindex(W.S₂.data)
    resize!(W.S₂.data[j][2],i)
    W.S₂.data[j][3] != nothing && resize!(W.S₂.data[j][3],i)
  end
  return nothing
end

function deleteat_stack!(W::NoiseProcess,i)
  for j in eachindex(W.S₂.data)
    deleteat!(W.S₂.data[j][2],i)
    W.S₂.data[j][3] != nothing && deleteat!(W.S₂.data[j][3],i)
  end
  return nothing
end

#=
function addat_stack!(W::NoiseProcess,i)
  for j in eachindex(W.S₂.data)
    resize!(W.S₂.data[j],i)
    W.S₂.data[j][3] != nothing && resize!(W.S₂.data[j][3],i)
  end
end
=#

iscontinuous(W::NoiseProcess) = W.continuous
