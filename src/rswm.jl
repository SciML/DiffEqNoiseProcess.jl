mutable struct RSWM{T}
  discard_length::T
  adaptivealg::Symbol
end

Base.@pure function RSWM(;
     discard_length=1e-15,
     adaptivealg::Symbol=:RSwM3)
     RSWM{typeof(discard_length)}(discard_length,adaptivealg)
end

adaptive_alg(rswm::RSWM{T}) where {T} = rswm.adaptivealg
