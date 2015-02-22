module TrustNLS

function maxStepLen(x::Vector, d::Vector, l::Vector, u::Vector)
  n = length(x)
  gamma = Inf
  for i = 1:n
    if d[i] > 0 && u[i] < Inf
      gamma = min(gamma, (u[i]-x[i])/d[i])
    elseif d[i] < 0 && l[i] > -Inf
      gamma = min(gamma, (l[i]-x[i])/d[i])
    end
  end
  return gamma
end

function solve(ncon::Int, h!::Function, J!::Function, x::Array{Cdouble,1},
    lower::Vector = [], upper::Vector = []; verbose::Bool = false)
  eps = 1e-4;
  nvar = length(x);
  theta = 0.9995;
  Delta = 1e3;
  beta1 = 0.1;
  beta2 = 0.25;
  kmax = 1000;

  if (length(lower) > 0 || length(upper) > 0)
    has_bounds = true
    length(lower) == 0 && (lower = -Inf*ones(nvar))
    length(upper) == 0 && (upper = Inf*ones(nvar))
  else
    has_bounds = false
  end

  hx = Array(Cdouble, ncon)
  hdn_plus = Array(Cdouble, ncon)
  A = Array(Cdouble, ncon, nvar)
  h!(x, hx);
  J!(x, A);

  verbose && println("x0 = ", x)
  verbose && println("h(x0) = ", hx)
  verbose && println("J(x0) = ", A)
  k = 0
  while k <= kmax
    verbose && println("|hx| = ", norm(hx))
    v = A'*hx;
    D = ones(nvar);
    if has_bounds
      for i = 1:nvar
        if v[i] < 0 && upper[i] < Inf
          D[i] = 1.0/sqrt(upper[i]-x[i]);
        elseif v[i] > 0 && lower[i] > -Inf
          D[i] = 1.0/sqrt(x[i]-lower[i]);
        end
      end
    end
    norm(v./D) < 1e-6 && break
    d = -v./(D.^2);

    alpha_cp = min(dot(v,d)/norm(A*d)^2, Delta/norm(D.*d));
    dcp = alpha_cp*d;
    has_bounds ? (gamma = maxStepLen(x, dcp, lower, upper)) : (gamma = Inf)
    gamma < 1.0 && (dcp = max(theta, 1-norm(dcp))*gamma*dcp)

    # This is not exactly what we want.
    # dn is an approximation to min |Ad+hx| s.t. |Dd| <= Delta
    dn = -A\hx;
    norm(D.*dn) > Delta && (dn = (Delta/norm(D.*dn))*dn)
    has_bounds ? (gamma = maxStepLen(x, dn, lower, upper)) : (gamma = Inf)
    gamma < 1.0 && (dn = max(theta, 1-norm(dn))*gamma*dn)

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
    h!(x+dn_plus, hdn_plus)
    rho_h = (nh2-norm(hdn_plus)^2)/(nh2-norm(hx+A*dn_plus)^2);

    if rho_h >= beta2
      x = x + dn_plus;
      verbose && println("x = ",x)
      Delta *= 2;
      h!(x, hx);
      J!(x, A);
    else
      Delta /= 4;
    end

    k += 1
    verbose && println("---")
  end

  return x, k, hx, A
end

end # module
