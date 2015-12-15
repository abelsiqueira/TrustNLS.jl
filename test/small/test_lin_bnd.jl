tests = 0
succ = 0
for nvar = 2:10
  for ncon = 2:10
    A = rand(ncon, nvar);
    while cond(A) > 1e4
      A = rand(ncon, nvar);
    end
    b = rand(ncon)
    l = rand(nvar)-2*rand(nvar)
    u = l+rand(nvar)+1
    function h!(x::Vector, h::Vector)
      return copy!(h, A*x-b)
    end

    function J!(x::Vector, J::Matrix)
      return copy!(J, A)
    end

    t = rand(nvar)
    x0 = t.*l+(1-t).*u;
    (x, k, hx, A) = TrustNLS.solve(ncon, h!, J!, x0, l, u)
    @test reduce(&, [l[i] <= x[i] <= u[i] for i = 1:nvar])
    tests += 1
    norm(bnd_scaling(A'*hx, x, l, u)) < 1e-6 && (succ += 1)
  end
end
println("Random Linear Bounded")
println("Success ratio: ", succ/tests)
@test succ/tests > 0.9
