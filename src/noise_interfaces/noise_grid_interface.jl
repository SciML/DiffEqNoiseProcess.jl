function save_noise!(W::NoiseGrid)

end

function linear_interpolant(Θ,dt,u0::Number,u1)
  (1-Θ)*u0 + Θ*u1
end

function linear_interpolant!(out,Θ,dt,u0,u1)
  Θm1 = (1-Θ); out .= Θm1.*u0 .+ Θ.*u1
end

function linear_interpolant(Θ,dt,u0::AbstractArray,u1)
  out = similar(u0)
  linear_interpolant!(out,Θ,dt,u0,u1)
  out
end

function interpolate!(W::NoiseGrid,t)
  ts,timeseries,timeseries2 = W.t,W.W,W.Z
  sign(W.dt)*t > sign(W.dt)*ts[end] && error("Solution interpolation cannot extrapolate past the final timepoint. Build a longer NoiseGrid to cover the integration.")
  sign(W.dt)*t < sign(W.dt)*ts[1] && error("Solution interpolation cannot extrapolate before the first timepoint. Build a longer NoiseGrid to cover the integration.")
  tdir = sign(ts[end]-ts[1])

  if t isa Union{Rational,Integer}
    @inbounds i = searchsortedfirst(ts,t,rev=tdir<0) # It's in the interval ts[i-1] to ts[i]
  else
    @inbounds i = searchsortedfirst(ts,t-tdir*10eps(typeof(t)),rev=tdir<0)
  end

  @inbounds if (t isa Union{Rational,Integer} && ts[i] == t) || (isapprox(t, ts[i]; atol = 100eps(typeof(t)), rtol = 100eps(t)))
    val1 = timeseries[i]
    timeseries2 != nothing ? val2 = timeseries2[i] : val2 = nothing
  elseif ts[i-1] == t # Can happen if it's the first value!
    val1 = timeseries[i-1]
    timeseries2 != nothing ? val2 = timeseries2[i-1] : val2 = nothing
  else
    dt = ts[i] - ts[i-1]
    Θ = (t-ts[i-1])/dt
    val1 = linear_interpolant(Θ,dt,timeseries[i-1],timeseries[i])
    timeseries2 != nothing ? val2 = linear_interpolant(Θ,dt,timeseries2[i-1],timeseries2[i]) : val2 = nothing
  end
  val1,val2
end

function interpolate!(out1,out2,W::NoiseGrid,t)
  ts,timeseries,timeseries2 = W.t,W.W,W.Z
  sign(W.dt)*t > sign(W.dt)*(ts[end]+10*sign(W.dt)*eps(typeof(t))) && error("Solution interpolation cannot extrapolate past the final timepoint. Build a longer NoiseGrid to cover the integration.")
  sign(W.dt)*t < sign(W.dt)*(ts[1]+10*sign(W.dt)*eps(typeof(t))) && error("Solution interpolation cannot extrapolate before the first timepoint. Build a longer NoiseGrid to cover the integration.")
  tdir = sign(ts[end]-ts[1])


  if t isa Union{Rational,Integer}
    @inbounds i = searchsortedfirst(ts,t,rev=tdir<0) # It's in the interval ts[i-1] to ts[i]
  else
    @inbounds i = searchsortedfirst(ts,t-tdir*10eps(typeof(t)),rev=tdir<0)
  end

  @inbounds if (t isa Union{Rational,Integer} && ts[i] == t) || (isapprox(t, ts[i]; atol = 100eps(typeof(t)), rtol = 100eps(t)))
    copyto!(out1,timeseries[i])
    timeseries2 != nothing && copyto!(out2,timeseries2[i])
  elseif ts[i-1] == t # Can happen if it's the first value!
    copyto!(out1,timeseries[i-1])
    timeseries2 != nothing && copyto!(out2,timeseries2[i-1])
  else
    dt = ts[i] - ts[i-1]
    Θ = (t-ts[i-1])/dt
    linear_interpolant!(out1,Θ,dt,timeseries[i-1],timeseries[i])
    timeseries2 != nothing && linear_interpolant!(out2,Θ,dt,timeseries2[i-1],timeseries2[i])
  end
end

function calculate_step!(W::NoiseGrid,dt,u,p)
  t = W.curt+dt
  if typeof(t) <: AbstractFloat && abs(t - W.t[end]) < 100eps(typeof(dt))
    t = W.t[end]
  end
  if isinplace(W)
    interpolate!(W.dW,W.dZ,W,t)
    W.dW .-= W.curW
    if W.Z != nothing
      W.dZ .-= W.curZ
    end
  else
    new_W, new_Z = W(t)
    W.dW = new_W - W.curW
    if W.Z != nothing
      W.dZ = new_Z - W.curZ
    end
  end
  W.dt = dt
end

function accept_step!(W::NoiseGrid,dt,u,p,setup_next=true)
  W.step_setup == false && error("Stepped past the defined domain for the NoiseGrid")

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
  if (W.dt isa Union{Rational,Integer})
    if sign(W.dt)*(W.curt + W.dt) > sign(W.dt)*W.t[end]
      setup_next = false
      W.step_setup = false
    end
  else
    if sign(W.dt)*(W.curt + W.dt) > sign(W.dt)*(W.t[end]+sign(W.dt)*10eps(typeof(dt)))
      setup_next = false
      W.step_setup = false
    end  
  end

  if setup_next
    calculate_step!(W,dt,u,p)
  end
  return nothing
end

function reject_step!(W::NoiseGrid,dtnew,u,p)
  calculate_step!(W,dtnew,u,p)
  return nothing
end

function setup_next_step!(W::NoiseGrid,u,p)
  calculate_step!(W,W.dt,u,p)
  return nothing
end
