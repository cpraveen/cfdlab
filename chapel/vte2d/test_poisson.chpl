// Solve
//   -Laplace(u) = f in Omega
//            u  = 0 on dOmega
use Math;
use StencilDist;
use Poisson;
use VTK;

config const n = 100,
             rtol = 1.0e-6,
             itmax = 1000;

const h = 1.0 / (n - 1);
const D = stencilDist.createDomain({1..n, 1..n}, fluff=(1,1));

// RHS function f in Poisson equation
proc rhsfun(x, y) { return sin(2*pi*x) * sin(2*pi*y); }

proc main()
{
   const x = [i in 1..n] (i-1)*h;
   var u, rhs : [D] real;
   forall (i,j) in D do
      rhs[i,j] = rhsfun(x[i], x[j]);

   const (res0,res,it) = poisson(u, rhs, h, rtol, itmax);
   writeln("Initial residual = ", res0);
   writeln("Final   residual = ", res);
   writeln("No. of iterations= ", it);

   write_vtk(x, x, 0.0, 0, "sol", u, "poisson.vtk");
}
