// Solves
//  -Laplace(u) = rhs
// u already has boundary values filled in which we do not change.
module Poisson
{
use Math;

private proc residual(u : [?D], rhs, h)
{
   const inner = D.expand(-1);
   var res = 0.0;
   forall (i,j) in inner with (+ reduce res)
   {
         const tmp = (u[i-1,j] - 2.0 * u[i,j] + u[i+1,j] 
                    + u[i,j-1] - 2.0 * u[i,j] + u[i,j+1]) / h**2
                    + rhs[i,j];
         res += tmp**2;
   }

   return sqrt(res/inner.size);
}

// red-black gauss-seidel
proc poisson(ref u : [?D], rhs, h, RTOL=1.0e-6, ITMAX=1000)
{
   const inner = D.expand(-1);
   const r = 2.0/(1.0 + pi * h);

   // Set initial guess to zero
   u = 0.0;

   var res0 = residual(u, rhs, h);
   var res  = res0;
   var it   = 0;

   while res > RTOL * res0 && it < ITMAX
   {
      forall i in inner.dim(0) do
      forall j in inner.dim(1) by 2 align i
      {
         const tmp = 0.25 * (u[i-1,j] + u[i+1,j] + u[i,j-1] + u[i,j+1] 
                             + h * h * rhs[i,j]);
         u[i,j] = (1.0 - r) * u[i,j] + r * tmp;
      }

      forall i in inner.dim(0) do
      forall j in inner.dim(1) by 2 align i+1
      {
         const tmp = 0.25 * (u[i-1,j] + u[i+1,j] + u[i,j-1] + u[i,j+1] 
                             + h * h * rhs[i,j]);
         u[i,j] = (1.0 - r) * u[i,j] + r * tmp;
      }

      res = residual(u, rhs, h);
      it += 1;
   }

   if res > RTOL * res0 && it == ITMAX
   {
      writeln("Error: no convergence in stream function");
      writeln("Error: increase RTOL and/or ITMAX");
      writef("res0, res, iter = %er %er %i\n", res0, res, it);
      exit();
   }

   return (res0, res, it);
}

}
