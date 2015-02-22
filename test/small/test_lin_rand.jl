for nvar = 5:5:50
  for ncon = 5:5:50
    A = rand(ncon, nvar);
    while rank(A) < min(ncon,nvar)
      A = rand(ncon, nvar);
    end
    b = rand(ncon)
    function h!(x::Vector, h::Vector)
      return copy!(h, A*x-b)
    end

    function J!(x::Vector, J::Matrix)
      return copy!(J, A)
    end

    (x, k, hx, A) = TrustNLS.solve(ncon, h!, J!, -rand(nvar))
    @test norm(A'*hx) < 1e-4
    @test k == 1
    y = A\b
    nvar < ncon && @test norm(x - y) < 1e-6
  end
end
