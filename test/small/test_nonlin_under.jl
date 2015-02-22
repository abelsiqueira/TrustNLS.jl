function h!(x::Vector, h::Vector)
  copy!(h, [x[2] - exp(x[1])])
end

function J!(x::Vector, J::Matrix)
  copy!(J, [-exp(x[1]) 1])
end

(x, k, hx, A) = TrustNLS.solve(1, h!, J!, [2.0;3.0])

@test norm(hx) < 1e-4
