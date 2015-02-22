using TrustNLS
using Base.Test

h(x) = [   x[1]^2+16*x[2]^2-17; 4*x[1]^2+4*x[2]^2-8]
J(x) = [ 2*x[1] 32*x[2]; 8*x[1]  8*x[2]]

sol = [1;1]

x = TrustNLS.solve(h, J, [2.0;3.0])
println("x = ", x)
