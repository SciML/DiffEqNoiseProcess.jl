@recipe function f(W::AbstractNoiseProcess)
  linewidth --> 3
  W.t,vecvec_to_mat(W.W)
end
