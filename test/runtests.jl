using TrustNLS
using Base.Test

function bnd_scaling(d::Vector, x::Vector, l::Vector, u::Vector)
  v = copy(d)
  for i = 1:length(v)
    if v[i] < 0 && u[i] < Inf
      v[i] *= sqrt(u[i] - x[i])
    elseif v[i] > 0 && l[i] > -Inf
      v[i] *= sqrt(x[i] - l[i])
    end
  end
  return v
end

include("small/test_nonlin_sq.jl")
include("small/test_nonlin_over.jl")
include("small/test_nonlin_under.jl")

include("small/test_lin_exact.jl")
include("small/test_lin_rand.jl")

include("small/test_bnd.jl")
