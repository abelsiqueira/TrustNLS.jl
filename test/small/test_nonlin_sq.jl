function h!(x::Vector, h::Vector)
  copy!(h, [x[1]^2+16*x[2]^2-17; 4*x[1]^2+4*x[2]^2-8]);
end

function J!(x::Vector, J::Matrix)
  copy!(J, [2*x[1] 32*x[2]; 8*x[1]  8*x[2]]);
end

sol = [1;1]

(x, k, hx, A) = TrustNLS.solve(2, h!, J!, [2.0;3.0])

@test norm(x-sol) < 1e-6
@test norm(A'*hx) < 1e-4
