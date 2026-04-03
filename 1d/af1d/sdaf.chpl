// Semi-discrete active flux
use StencilDist;
use Math;
use IO;
use PDE;

const ark : [1..3] real = (0.0, 3.0/4.0, 1.0/3.0);
const brk : [1..3] real = (1.0, 1.0/4.0, 2.0/3.0);

config const n    = 100,
             Tf   = 1.0,
             cfl  = 0.25,
             diff = 1;

const xmin = 0.0;
const xmax = 1.0;
const dx = (xmax - xmin)/n;
const xc = [i in 1..n] xmin + (i-0.5)*dx; // cell centers
const xv = [i in 1..n] xmin + (i-1.0)*dx; // vertices

const D = stencilDist.createDomain({1..n}, fluff=(1,), periodic=true);

//-----------------------------------------------------------------------------
proc initial_condition(x)
{
   return sin(2*pi*x);
}

//-----------------------------------------------------------------------------
proc compute_dt(u : real) : real
{
   const a = jacobian(u);
   return cfl * dx / (abs(a) + 1.0e-16);
}

//-----------------------------------------------------------------------------
// vertex = [i]
//         [i-1]--(i-1)--[i]--(i)--[i+1]
proc rhsv1(uc, uv, ref Rv)
{
   forall i in D
   {
      const a = jacobian(uv[i]);
      var ux : real;
      if a > 0.0
      {
         // backward difference
         ux = (2.0 * uv[i-1] - 6.0 * uc[i-1] + 4.0 * uv[i]) / dx;
      }
      else
      {
         // forward difference
         ux = (-4.0 * uv[i] + 6.0 * uc[i] - 2.0 * uv[i+1]) / dx ;
      }
      Rv[i] = a * ux;
   }
}

//-----------------------------------------------------------------------------
proc rhsc(uv, ref Rc)
{
   forall i in D
   {
      const fl = flux(uv[i]);
      const fr = flux(uv[i+1]);
      Rc[i] = (fr - fl) / dx;
   }
}

//-----------------------------------------------------------------------------
proc output(x, u, filename)
{
   var fw = open(filename, ioMode.cw).writer(locking=false);
   for i in D
   {
      fw.writef("%14.8er %14.8er\n", x[i], u[i]);
   }
   fw.close();
}

//-----------------------------------------------------------------------------
proc main()
{
   var uc,  uv  : [D] real;
   var uc0, uv0 : [D] real;
   var Rc,  Rv  : [D] real;

   // Set initial condition
   forall i in 1..n
   {
      uv[i] = initial_condition(xv[i]);
   }
   uv.updateFluff();

   forall i in 1..n
   {
      const w = initial_condition(xc[i]);
      uc[i] = (uv[i] + 4.0 * w + uv[i+1]) / 6.0;
   }
   uc.updateFluff();

   output(xc, uc, "uc0.txt");
   output(xv, uv, "uv0.txt");

   var t = 0.0;
   var dt = 0.0;
   while t < Tf
   {
      dt = min reduce compute_dt(uc);
      if t + dt > Tf then
         dt = Tf - t;

      uc0 = uc;
      uv0 = uv;

      for rk in 1..3
      {
         rhsv1(uc, uv, Rv);
         rhsc(uv, Rc);

         forall i in D
         {
            uc[i] = ark[rk] * uc0[i] + brk[rk] * (uc[i] - dt * Rc[i]);
            uv[i] = ark[rk] * uv0[i] + brk[rk] * (uv[i] - dt  * Rv[i]);
         }
         uc.updateFluff();
         uv.updateFluff();
      }
      t = t + dt;
      writeln("t, dt = ", t, dt);
   }

   output(xc, uc, "uc.txt");
   output(xv, uv, "uv.txt");
}
