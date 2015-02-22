for nvar = 5:5:50
  y = ones(nvar)
  for ncon = 5:5:50
    A = rand(ncon, nvar);
    b = A*y
    function h!(x::Vector, h::Vector)
      return copy!(h, A*x-b)
    end

    function J!(x::Vector, J::Matrix)
      return copy!(J, A)
    end

    (x, k, hx, A) = TrustNLS.solve(ncon, h!, J!, -rand(nvar))
    @test norm(A'*hx) < 1e-4
    @test k == 1
    @test norm(hx) < 1e-4
    nvar < ncon && @test norm(x - y) < 1e-6
  end
end
