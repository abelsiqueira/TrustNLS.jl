using CUTEst
using TrustNLS

if length(ARGS) == 0
  println("Need CUTEst problem as argument")
  exit(1)
end

curdir = pwd();
tmpdir = mktempdir();
cd(tmpdir)

nlp = CUTEstModel(convert(ASCIIString,ARGS[1]))

if nlp.meta.ncon == 0
  println("Problem is unconstrained, EXITING")
  exit(1)
end

function h!(x::Vector, h::Vector)
  st = Int32[0]
  nnzj = Int32[0]
  CUTEst.ccfsg(st, Int32[nlp.meta.nvar], Int32[nlp.meta.ncon], x, h, nnzj,
      Int32[nlp.meta.nnzj], [0.0], Int32[0], Int32[0], Int32[0], nlp.libname)
end

function J!(x::Vector, J::Matrix)
  st = Int32[0]
  nnzj = Int32[0]
  h = Array(Cdouble, nlp.meta.ncon)
  jval = Array(Cdouble, nlp.meta.nnzj)
  jvar = Array(Cint, nlp.meta.nnzj)
  jfun = Array(Cint, nlp.meta.nnzj)
  CUTEst.ccfsg(st, Int32[nlp.meta.nvar], Int32[nlp.meta.ncon], x, h, nnzj,
      Int32[nlp.meta.nnzj], jval, jvar, jfun, Int32[1], nlp.libname)
  for i = 1:nlp.meta.ncon
    for j = 1:nlp.meta.nvar
      J[i,j] = 0.0
    end
  end
  for k = 1:nnzj[1]
    J[jfun[k],jvar[k]] = jval[k]
  end
end

@time (x,k,hx,Jx) = TrustNLS.solve(nlp.meta.ncon, h!, J!, nlp.meta.x0)

println("|hx| = ", norm(hx))
println("k = ", k)

cd(curdir)
