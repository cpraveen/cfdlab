use Math;
use StencilDist;
use Poisson;
use VTK;

config const n = 128,
             Re = 100.0,
             Tf = 100.0,
             wtol = 1.0e-4,
             ptol = 1.0e-4,
             witmax = 50000,
             pitmax = 1000;

const nu = 1.0/Re;
const h = 1.0/(n - 1);
const D = stencilDist.createDomain({1..n, 1..n}, fluff=(1,1));
const inner = D.expand(-1);

//-----------------------------------------------------------------------------
// Compute time step based on Fourier stability
//-----------------------------------------------------------------------------
proc time_step(u, v)
{
   var dt = 0.25 * h**2 / nu;
   forall ij in D with (min reduce dt)
   {
      dt = min(dt, 2.0 * nu /(u[ij]**2 + v[ij]**2 + 1.0e-12));
   }

   return dt;
}

//-----------------------------------------------------------------------------
// Compute velocit from stream function
//-----------------------------------------------------------------------------
proc compute_velocity(psi, ref u, ref v)
{
   forall (i,j) in inner
   {
      u[i,j] =  (psi[i,j+1] - psi[i,j-1]) / (2*h);
      v[i,j] = -(psi[i+1,j] - psi[i-1,j]) / (2*h);
   }

   u.updateFluff();
   v.updateFluff();
}

//-----------------------------------------------------------------------------
// Fill boundary values of omega
//-----------------------------------------------------------------------------
proc boundary(psi, u, v, ref omega)
{
   const ih2 = 2.0 / h**2;

   forall i in 2..n-1
   {
      // bottom boundary
      omega[i,1] = ih2 * (psi[i,1] - psi[i,2] + h * u[i,1]);

      // top boundary
      omega[i,n] = ih2 * (psi[i,n] - psi[i,n-1] - h * u[i,n]);
   }

   forall j in 2..n-1
   {
      // left boundary
      omega[1,j] = ih2 * (psi[1,j] - psi[2,j] - h * v[1,j]);

      // right boundary
      omega[n,j] = ih2 * (psi[n,j] - psi[n-1,j] + h * v[n,j]);
   }
}

//-----------------------------------------------------------------------------
proc update_vort(dt, psi, u, v, ref omega)
{
   const w = omega; // make copy
   var wres = 0.0;

   forall (i,j) in inner with (+ reduce wres)
   {
      const wx = (w[i+1,j] - w[i-1,j]) / (2*h);
      const wy = (w[i,j+1] - w[i,j-1]) / (2*h);
      const wxx = (w[i-1,j] - 2 * w[i,j] + w[i+1,j]) / h**2;
      const wyy = (w[i,j-1] - 2 * w[i,j] + w[i,j+1]) / h**2;
      omega[i,j] = w[i,j] + dt*( - u[i,j] * wx - v[i,j] * wy 
                                 + nu * (wxx + wyy) );
      wres += ((omega[i,j] - w[i,j])/dt)**2;
   }

   omega.updateFluff();

   return sqrt(wres/inner.size);
}

//-----------------------------------------------------------------------------
proc main()
{
   const x : [1..n] real = [i in 1..n] (i-1)*h;

   var psi, omega, u, v : [D] real;

   // Velocity at top lid
   u[..,n] = 1.0;

   var t = 0.0, it = 0, wres = 1.0e20;
   while t < Tf && it < witmax && wres > wtol
   {
      const (res0,res,pit) = poisson(psi, omega, h, ptol, pitmax);
      writef("stream: res0,res,it      = %10.3er %10.3er %4i\n",res0,res,pit);
      compute_velocity(psi, u, v);
      boundary(psi, u, v, omega);
      const dt = time_step(u, v);
      wres = update_vort(dt, psi, u, v, omega);
      t += dt; it += 1;
      const omin = min reduce omega;
      const omax = max reduce omega;
      writef("vort  : it,t,res,min,max = %4i %10.3er %10.3er %10.3er %10.3er\n",
              it,t,wres,omin,omax);
   }

   const (res0,res,pit) = poisson(psi, omega, h);
   compute_velocity(psi, u, v);
   boundary(psi, u, v, omega);
   const fname = "sol.vtk";
   write_vtk(x, x, t, it, ("psi","omega","u","v"), (psi,omega,u,v), fname,
             "Re = " + Re:int:string);
}
