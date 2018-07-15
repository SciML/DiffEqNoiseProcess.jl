mutable struct RSWM{adaptivealg,T}
  discard_length::T
end

Base.@pure function RSWM(;
     discard_length=1e-15,
     adaptivealg::Symbol=:RSwM3)
     RSWM{adaptivealg,typeof(discard_length)}(discard_length)
end

adaptive_alg(rswm::RSWM{adaptivealg,T}) where {adaptivealg,T} = adaptivealg
