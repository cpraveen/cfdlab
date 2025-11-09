// Solve -u'' = f

module Poisson
{

use Math;

//------------------------------------------------------------------------------
// Computes residual norm
//------------------------------------------------------------------------------
private
proc residual(v : [?D] real,
              f : [D] real,
              h : real,
              ref r : [D] real)
{
   const n = D.dim(0).size;
   const inner = D.expand(-1);

   // Due to dirichlet bc
   r[1] = 0.0;
   r[n] = 0.0;

   var rnorm = 0.0;
   forall i in inner with (+ reduce rnorm)
   {
      r[i] = f[i] + (v[i-1] - 2.0 * v[i] + v[i+1]) / h**2;
      rnorm += r[i]**2;
   }

   return sqrt(rnorm/inner.size);
}

//------------------------------------------------------------------------------
// Weighted jacobi iterations
//------------------------------------------------------------------------------
private
proc wjacobi(ref v : [?D], f, h, niter)
{
   const inner = D.expand(-1);
   const w = 2.0/3.0;

   var vold = v;

   for it in 1..niter
   {      
      vold <=> v;
      forall i in inner
      {
         v[i] = 0.5 * (vold[i-1] + vold[i+1] + h**2 * f[i]);
         v[i] = (1.0-w) * vold[i] + w * v[i];
      }
   }

}

//------------------------------------------------------------------------------
// vh --> v2h
//------------------------------------------------------------------------------
private
proc restrictfw(vh : [?Dh], 
                ref v2h : [?D2h])
{
   // odd indices in D2h: inject
   forall i in D2h by 2 align 1
   {
      v2h[i] = vh[2*i-1];
   }

   // even indices in D2h: average
   forall i in D2h by 2 align 2
   {
      const j = 2 * i - 1;
      v2h[i] = 0.25 * (vh[j-1] + 2 * vh[j] + vh[j+1]);
   }
}

//------------------------------------------------------------------------------
// v2h --> vh
//------------------------------------------------------------------------------
private
proc  prolongate(v2h : [?D2h], 
                 ref vh : [?Dh])
{
   // Injection for odd indices in Dh: 1,3,5,...
   forall i in D2h do
      vh[2*i-1] = v2h[i];

   // Average for even indices in Dh: 2,4,6,...
   forall i in Dh by 2 align 2 do
      vh[i] = 0.5 * (vh[i-1] + vh[i+1]);
}

//------------------------------------------------------------------------------
// See Briggs et al., page 40
//------------------------------------------------------------------------------
private
proc vcycle(ref vh : [?Dh], 
            fh : [Dh], 
            h : real,
            nsmooth : int,
            levels : int,
            Lh : int) : real
{
   wjacobi(vh, fh, h, nsmooth);
   var rh : [Dh] real;

   if Lh != levels // not on coarsest grid
   {
      residual(vh, fh, h, rh);

      const nh = Dh.dim(0).size;
      const n2h = (nh + 1) / 2 : int;
      const D2h = {1..n2h};
      var f2h : [D2h] real;
      restrictfw(rh, f2h);

      var e2h : [D2h] real;
      const r2hnorm = vcycle(e2h, f2h, 2*h, nsmooth, levels, Lh+1);

      var eh : [Dh] real;
      prolongate(e2h, eh);
      vh += eh;
   }
   
   wjacobi(vh, fh, h, nsmooth);
   const rhnorm = residual(vh, fh, h, rh);
   return rhnorm;
}

//------------------------------------------------------------------------------
// v-cycle multigrid
// v should have boundary values filled in.
//------------------------------------------------------------------------------
proc multigrid(ref v : [?D], f, h, rtol, niter, nsmooth, levels)
{
   const inner = D.expand(-1);

   var fnorm = 0.0;
   forall i in inner with (+ reduce fnorm) do
      fnorm += f[i]**2;
   fnorm = sqrt(fnorm / inner.size);

   var r : [D] real;
   var rnorm = residual(v, f, h, r);

   var it = 0;
   while rnorm > rtol * fnorm && it < niter
   {
      const rnorm_new = vcycle(v, f, h, nsmooth, levels, 1);
      const conv = rnorm_new / rnorm;
      rnorm = rnorm_new;
      it += 1;
      writef("it,rnorm,conv = %4i %12.4er %12.4er\n", it, rnorm, conv);
   }

   if rnorm > rtol * fnorm
   {
      writeln("Specified tolerance not achieved !!!");
      writeln("Increase niter and rerun");
   }

}

//------------------------------------------------------------------------------
// Red-black SOR
// v should have boundary values filled in.
//------------------------------------------------------------------------------
proc sor(ref v : [?D], f, h, rtol, niter)
{
   const omg = 2.0 / (1.0 + sin(pi*h));
   const inner = D.expand(-1);

   var fnorm = 0.0;
   forall i in inner with (+ reduce fnorm) do
      fnorm += f[i]**2;
   fnorm = sqrt(fnorm / inner.size);

   var r : [D] real;
   var rnorm = residual(v, f, h, r);

   var it = 0;
   while rnorm > rtol * fnorm && it < niter
   {
      forall i in inner by 2 align 1
      {
         const tmp = 0.5 * (v[i-1] + v[i+1] + h**2 * f[i]);
         v[i] = (1.0 - omg) * v[i] + omg * tmp;
      }

      forall i in inner by 2 align 2
      {
         const tmp = 0.5 * (v[i-1] + v[i+1] + h**2 * f[i]);
         v[i] = (1.0 - omg) * v[i] + omg * tmp;
      }

      const rnorm_new = residual(v, f, h, r);
      const conv = rnorm_new / rnorm;
      rnorm = rnorm_new;
      it += 1;
      writef("it,rnorm,conv = %4i %12.4er %12.4er\n", it, rnorm, conv);
   }

   if rnorm > rtol * fnorm
   {
      writeln("Specified tolerance not achieved !!!");
      writeln("Increase niter and rerun");
   }

}

}
