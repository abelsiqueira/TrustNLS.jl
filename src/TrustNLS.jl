module TrustNLS

function solve(h::Function, J::Function, x::Array{Cdouble,1},
    verbose::Bool = false)
  eps = 1e-4;
  n = length(x);
  theta = 0.9995;
  Delta = 1e3;
  beta1 = 0.1;
  beta2 = 0.25;
  kmax = 1000;

  verbose && println("x0 = ", x)
  hx = h(x);
  A = J(x);
  k = 0
  while norm(A*hx) > eps
    verbose && println("|hx| = ", norm(hx))
    k += 1
    if k > kmax
      break
    end
    v = A'*hx;
    D = ones(n);
    for i = 1:n
      if v[i] < 0
        D[i] = 1.0/sqrt(Delta-x[i]);
      elseif v[i] > 0
        D[i] = 1.0/sqrt(x[i]+Delta);
      end
    end
    d = -v./(D.^2);

    alpha_cp = min(dot(v,d)/norm(A*d)^2, Delta/norm(D.*d));
    dcp = alpha_cp*d;
    outside = false;
    gamma = 1.0;
    for i = 1:n
      if dcp[i] > 0
        gamma = min(gamma, Delta/dcp[i]);
        outside = true;
      elseif dcp[i] < 0
        gamma = min(gamma, -Delta/dcp[i]);
        outside = true;
      end
    end
    if outside
      dcp = max(theta, 1-norm(dcp))*gamma*dcp;
    end

    dn = -A\hx;
    if norm(D.*dn) > Delta
      dn = (Delta/norm(D.*dn))*dn
    end
    outside = false;
    gamma = 1.0
    for i = 1:n
      if dn[i] > 0
        gamma = min(gamma, Delta/dn[i]);
        outside = true;
      elseif dn[i] > 0
        gamma = min(gamma, -Delta/dn[i]);
        outside = true;
      end
    end
    if outside
      dn = max(theta, 1-norm(dn))*gamma*dn;
    end

    verbose && println("dcp = ",dcp)
    verbose && println("dn = ",dn)
    nh2 = norm(hx)^2;
    nhdcp2 = norm(hx + A*dcp)^2;
    rho_C(delta) = (nh2 - norm(hx+A*delta)^2)/(nh2 - nhdcp2);
    t = 1.0;
    while rho_C(t*dn+(1-t)*dcp) >= beta1
      t = 0.9*t;
      if (t < 1e-6)
        verbose && println("Warning")
        break
      end
    end
    dn_plus = t*dn + (1-t)*dcp;
    rho_h = (nh2-norm(h(x+dn_plus))^2)/(nh2-norm(hx+A*dn_plus)^2);

    if rho_h >= beta2
      x = x + dn_plus;
      verbose && println("x = ",x)
      Delta *= 2;
      hx = h(x);
      A = J(x);
    else
      Delta /= 4;
    end

    verbose && println("---")
  end

  return x
end

end # module
