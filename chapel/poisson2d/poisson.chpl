module Poisson
{

use Math;

//------------------------------------------------------------------------------
record MGData
{
   const D : domain(2);
   const dx, dy : real;
   var v, f, r : [D] real;
}

//------------------------------------------------------------------------------
private proc residual(v : [?D], f, ref r, dx, dy)
{
   const inner = D.expand(-1);
   const rdx2 = 1.0/dx**2;
   const rdy2 = 1.0/dy**2;
   const nx = D.dim(0).last,
         ny = D.dim(1).last;

   r[0 , ..] = 0.0;
   r[nx, ..] = 0.0;
   r[.., 0 ] = 0.0;
   r[.., ny] = 0.0;

   var rnorm = 0.0;
   forall (i,j) in inner with (+ reduce rnorm)
   {
      r[i,j] =   rdx2 * (v[i-1,j] - 2 * v[i,j] + v[i+1,j]) 
               + rdy2 * (v[i,j-1] - 2 * v[i,j] + v[i,j+1]) 
               + f[i,j];
      rnorm += r[i,j]**2;
   }

   return sqrt(rnorm/inner.size);
}

//------------------------------------------------------------------------------
private proc residual(ref data)
{
   const dx = data.dx;
   const dy = data.dy;

   ref v = data.v;
   ref f = data.f;
   ref r = data.r;

   return residual(v, f, r, dx, dy);
}

//------------------------------------------------------------------------------
// Weighted Jacobi
//------------------------------------------------------------------------------
private proc wjacobi(ref data, niter)
{
   const omg = 4.0/5.0; // optimal weight for 2d poisson
   const dx    = data.dx, 
         dy    = data.dy,
         rdx2  = 1.0/dx**2, 
         rdy2  = 1.0/dy**2,
         rf    = 1.0/(2.0/dx**2 + 2.0/dy**2),
         inner = data.D.expand(-1);

   ref v = data.v;
   ref f = data.f;
   var vold = v;

   for n in 1..niter
   {
      vold <=> v;
      forall (i,j) in inner
      {
         v[i,j] = rf * (  rdx2*(vold[i-1,j] + vold[i+1,j]) 
                        + rdy2*(vold[i,j-1] + vold[i,j+1]) 
                        + f[i,j]);
         v[i,j] = (1.0-omg) * vold[i,j] + omg * v[i,j];
      }
   }
}

//------------------------------------------------------------------------------
// vh --> v2h
// Full weighting, Briggs, page 36
// Boundary values of v2h are assumed to be already zero.
//------------------------------------------------------------------------------
private proc restrict(vh : [?Dh],
                      ref v2h : [?D2h])
{
   const I2h = D2h.expand(-1);

   forall (i2,j2) in I2h
   {
      const i = 2 * i2;
      const j = 2 * j2;
      v2h[i2,j2] = ( vh[i-1,j-1] + vh[i-1,j+1] + vh[i+1,j-1] + vh[i+1,j+1]
                   + 2.0 * (vh[i,j-1] + vh[i,j+1] + vh[i-1,j] + vh[i+1,j])
                   + 4.0 * vh[i,j]) / 16.0;
   }
}

//------------------------------------------------------------------------------
// v2h --> vh
// Briggs, page 35
// Boundary values of vh are assumed to be already zero.
//------------------------------------------------------------------------------
private proc prolong(v2h : [?D2h],
                     ref vh : [?Dh])
{
   const nx  = D2h.dim(0).last;
   const ny  = D2h.dim(1).last;
   const S2h = {0..nx-1, 0..ny-1};

   forall (i,j) in S2h
   {
      vh[2*i  , 2*j  ] = v2h[i,j];
      vh[2*i+1, 2*j  ] = 0.5 * (v2h[i,j] + v2h[i+1,j]);
      vh[2*i  , 2*j+1] = 0.5 * (v2h[i,j] + v2h[i,j+1]);
      vh[2*i+1, 2*j+1] = 0.25 * (v2h[i,j]   + v2h[i+1,j] + 
                                 v2h[i,j+1] + v2h[i+1,j+1]);
   }
}

//------------------------------------------------------------------------------
private proc vcycle(ref data, nsmooth, L) : real
{
   wjacobi(data[L], nsmooth);

   if L != data.size
   {
      residual(data[L]);
      ref vh  = data[L].v, 
          rh  = data[L].r, 
          f2h = data[L+1].f, 
          v2h = data[L+1].v;
      restrict(rh, f2h);

      v2h = 0.0;
      const r2hnorm = vcycle(data, nsmooth, L+1);
      prolong(v2h, rh);
      vh += rh;
   }

   wjacobi(data[L], nsmooth);
   return residual(data[L]);
}

//------------------------------------------------------------------------------
// Domain D must be of the form {0..nx, 0..ny}
//------------------------------------------------------------------------------
proc multigrid(ref v : [?D], f, dx, dy, rtol, niter, levels, nsmooth)
{
   assert(D.dim(0).first == 0 &&
          D.dim(1).first == 0,
          "Starting index must be 0");

   const nx = D.dim(0).last;
   const ny = D.dim(1).last;

   assert(nx % 2**levels == 0 &&
          ny % 2**levels == 0, "nx, ny must be divisible by 2^levels");

   // Build multigrid hierarchy
   var NX, NY : [1..levels] int;
   var DX, DY : [1..levels] real;
   for l in 1..levels
   {
      NX[l] = nx / 2**(l-1);
      NY[l] = ny / 2**(l-1);
      DX[l] = 2**(l-1) * dx;
      DY[l] = 2**(l-1) * dy;
   }

   var data = [i in 1..levels] new MGData({0..NX[i],0..NY[i]}, 
                                          DX[i], DY[i]);
   data[1].v = v;
   data[1].f = f;

   const inner = D.expand(-1);
   var fnorm = + reduce [ij in inner] f[ij]**2;
   fnorm = sqrt(fnorm/inner.size);

   var rnorm = residual(data[1]);
   var it = 0;
   while rnorm > rtol * fnorm && it < niter
   {
      const rnorm_new = vcycle(data, nsmooth, 1);
      const rate = rnorm_new / rnorm;
      rnorm = rnorm_new;
      it += 1;
      writef("it,rnorm,rate = %5i %12.4er %12.4er\n",it,rnorm,rate);
   }

   v = data[1].v;
}

//------------------------------------------------------------------------------
// Red-black Gauss-Seidel
//------------------------------------------------------------------------------
proc sor(ref v : [?D], f, dx, dy, rtol, niter)
{
   const rdx2  = 1.0/dx**2, 
         rdy2  = 1.0/dy**2,
         rf    = 1.0/(2.0/dx**2 + 2.0/dy**2),
         h     = min(dx,dy),
         omg   = 2.0/(1.0 + sin(pi*h)),
         inner = D.expand(-1);

   var fnorm = + reduce [ij in inner] f[ij]**2;
   fnorm = sqrt(fnorm/inner.size);

   var r : [D] real;
   var rnorm = residual(v, f, r, dx, dy);
   var it = 0;
   while rnorm > rtol * fnorm && it < niter
   {
      forall (i,j) in inner do
      if (i+j)%2 == 0
      {
         const tmp = rf * (  rdx2*(v[i-1,j] + v[i+1,j]) 
                           + rdy2*(v[i,j-1] + v[i,j+1]) 
                           + f[i,j]);
         v[i,j] = (1.0-omg) * v[i,j] + omg * tmp;
      }

      forall (i,j) in inner do
      if (i+j)%2 == 1
      {
         const tmp = rf * (  rdx2*(v[i-1,j] + v[i+1,j]) 
                           + rdy2*(v[i,j-1] + v[i,j+1]) 
                           + f[i,j]);
         v[i,j] = (1.0-omg) * v[i,j] + omg * tmp;
      }

      const rnorm_new = residual(v, f, r, dx, dy);
      const rate = rnorm_new / rnorm;
      rnorm = rnorm_new;
      it += 1;
      writef("it,rnorm,rate = %5i %12.4er %12.4er\n",it,rnorm,rate);
   }

}

}
