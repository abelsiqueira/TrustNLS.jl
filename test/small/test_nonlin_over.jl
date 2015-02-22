function h!(x::Vector, h::Vector)
  return copy!(h, [x[1]^2-x[2]; x[1]*x[2]-4; x[1]^2+x[2]^2+1])
end

function J!(x::Vector, J::Matrix)
  return copy!(J, [2*x[1] -1; x[2] x[1]; 2*x[1] 2*x[2]])
end

(x, k, hx, A) = TrustNLS.solve(3, h!, J!, [2.0;3.0])

@test norm(A'*hx) < 1e-4
