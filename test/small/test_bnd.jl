function h!(x::Vector, h::Vector)
  return copy!(h, [x[1]-1; 10*(x[2]-x[1]^2)])
end

function J!(x::Vector, J::Matrix)
  return copy!(J, [1 0; -20*x[1] 10])
end

l = [-Inf;0.5]
u = [0.5;Inf]

(x, k, hx, A) = TrustNLS.solve(2, h!, J!, [0.0;0.5], l, u)
@test norm(bnd_scaling(A'*hx, x, l, u)) < 1e-4
@test reduce(&, [l[i] <= x[i] <= u[i] for i = 1:2])
