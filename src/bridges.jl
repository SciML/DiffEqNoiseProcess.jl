function BrownianBridge(t0,tend,W0,Wend,Z0=nothing,Zend=nothing;kwargs...)
  W = WienerProcess(t0,W0,Z0;kwargs...)
  h = tend-t0
  Wh = Wend-W0
  if Z0 != nothing
    Zh = Zend - Z0
  else
    Zh = nothing
  end
  push!(W.S₁,(h,Wh,Zh))
  W
end

function BrownianBridge!(t0,tend,W0,Wh,Z0=nothing,Zh=nothing;kwargs...)
  W = WienerProcess!(t0,W0,Z0;kwargs...)
  h = tend-t0
  Wh .-= W0
  if Z0 != nothing
    Zh .-= Z0
  else
    Zh = nothing
  end
  push!(W.S₁,(h,Wh,Zh))
  W
end

function GeometricBrownianBridge(μ,σ,t0,tend,W0,Wend,Z0=nothing,Zend=nothing;kwargs...)
  W = GeometricBrownianMotionProcess(μ,σ,t0,W0,Z0;kwargs...)
  h = tend-t0
  Wh = Wend-W0
  if Z0 != nothing
    Zh = Zend - Z0
  else
    Zh = nothing
  end
  push!(W.S₁,(h,Wh,Zh))
  W
end

function GeometricBrownianBridge!(μ,σ,t0,tend,W0,Wh,Z0=nothing,Zh=nothing;kwargs...)
  W = GeometricBrownianMotionProcess!(μ,σ,t0,W0,Z0;kwargs...)
  h = tend-t0
  Wh .-= W0
  if Z0 != nothing
    Zh .-= Z0
  else
    Zh = nothing
  end
  push!(W.S₁,(h,Wh,Zh))
  W
end
