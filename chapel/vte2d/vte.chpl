use Math;
use Poisson;
use VTK;

config const n = 100,
             Re = 10.0,
             Tf = 100.0,
             wtol = 1.0e-4,
             ptol = 1.0e-4,
             witmax = 1000,
             pitmax = 1000;

const nu = 1.0/Re;
const h = 1.0/(n - 1);
const D = {1..n, 1..n};
const inner = D.expand(-1);

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
}

//-----------------------------------------------------------------------------
// Fill boundary values of omega
//-----------------------------------------------------------------------------
proc boundary(psi, u, v, ref omega)
{
   forall i in 1..n
   {
      var u1, u2, u3 : real;

      // bottom boundary
      u1 = u[i,1];
      u2 = (psi[i,2] - psi[i,1]) / h;
      u3 = (psi[i,3] - psi[i,2]) / h;
      omega[i,1] = -(- 8 * u1 + 9 * u2 - u3) / (3*h);

      // top boundary
      u1 = u[i,n];
      u2 = (psi[i,n] - psi[i,n-1]) / h;
      u3 = (psi[i,n-1] - psi[i,n-2]) / h;
      omega[i,n] = -( 8 * u1 - 9 * u2 + u3) / (3*h);
   }

   forall j in 1..n
   {
      var v1, v2, v3 : real;

      // left boundary
      v1 = v[1,j];
      v2 = -(psi[2,j] - psi[1,j]) / h;
      v3 = -(psi[3,j] - psi[2,j]) / h;
      omega[1,j] = (- 8 * v1 + 9 * v2 - v3) / (3*h);

      // right boundary
      v1 = v[n,j];
      v2 = -(psi[n,j] - psi[n-1,j]) / h;
      v3 = -(psi[n-1,j] - psi[n-2,j]) / h;
      omega[n,j] = ( 8 * v1 - 9 * v2 + v3) / (3*h);
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
      wres += (omega[i,j] - w[i,j])**2;
   }

   return sqrt(wres/inner.size);
}

//-----------------------------------------------------------------------------
proc main()
{
   const x : [1..n] real = [i in 1..n] (i-1)*h;

   var psi, omega, u, v : [D] real;

   // Velocity at top lid
   u[..,n] = 1.0;

   const dt = 0.25 * Re * h**2;

   var t = 0.0;
   var it = 0;
   var wres = 1.0e20;
   while t < Tf && it < witmax && wres > wtol
   {
      const (res0,res,pit) = poisson(psi, omega, h, ptol, pitmax);
      writef("stream: res0,res,it      = %10.3er %10.3er %4i\n",res0,res,it);
      compute_velocity(psi, u, v);
      boundary(psi, u, v, omega);
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
   write_vtk(x, x, t, it, ["psi"], psi, fname);
   write_vtk(["omega"], omega, fname);
}
