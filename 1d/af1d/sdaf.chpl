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
             diff = 1,
             ic   = 1;

const xmin = 0.0;
const xmax = 1.0;
const dx = (xmax - xmin)/n;
const xc = [i in 1..n] xmin + (i-0.5)*dx; // cell centers
const xv = [i in 1..n] xmin + (i-1.0)*dx; // vertices

const D = stencilDist.createDomain({1..n}, fluff=(1,), periodic=true);

//-----------------------------------------------------------------------------
proc ic1(x)
{
   return sin(2*pi*x);
}

//-----------------------------------------------------------------------------
// Square hat
//-----------------------------------------------------------------------------
proc ic2(x)
{
   if abs(x-0.5) < 0.25
   {
      return 1.0;
   }
   else
   {
      return -0.5;
   }
}

//-----------------------------------------------------------------------------
// Square hat, continuous, spread over one dx
//-----------------------------------------------------------------------------
proc ic3(x)
{
   const w = 0.25, xm = 0.5;
   const u1 = 1.0, u2 = -0.5;

   if abs(x-xm) < w
   {
      return u1;
   }
   else if x > xm + w + dx || x < xm - (w + dx)
   {
      return u2;
   }
   else if x >= xm + w // in [xm+w, xm+w+dx]
   {
      const f = (x - (xm + w))/dx;
      return (1-f)*u1 + f*u2;
   }
   else // in [xm-w-dx,xm-w]
   {
      const f = (x - (xm - w - dx))/dx;
      return (1-f)*u2 + f*u1;
   }
}

//-----------------------------------------------------------------------------
proc initial_condition(x)
{
   if ic == 1 then
      return ic1(x);
   else if ic == 2 then
      return ic2(x);
   else
      return ic3(x);
}

//-----------------------------------------------------------------------------
proc compute_dt(u : real) : real
{
   const J = jacobian(u);
   return cfl * dx / (abs(J) + 1.0e-16);
}

//-----------------------------------------------------------------------------
// RHS at vertices; uses one sided derivative
// vertex = [i]
//         [i-1]--(i-1)--[i]--(i)--[i+1]
//-----------------------------------------------------------------------------
proc rhsv1(uc, uv, ref Rv)
{
   forall i in D
   {
      const J = jacobian(uv[i]);
      var ux : real;
      if J > 0.0
      {
         // backward difference
         ux = (2.0 * uv[i-1] - 6.0 * uc[i-1] + 4.0 * uv[i]) / dx;
      }
      else
      {
         // forward difference
         ux = (-4.0 * uv[i] + 6.0 * uc[i] - 2.0 * uv[i+1]) / dx ;
      }
      Rv[i] = J * ux;
   }
}

//-----------------------------------------------------------------------------
// RHS at vertices; uses stencil with one downwind bias
// vertex = [i]
//         [i-1]--(i-1)--[i]--(i)--[i+1]
//-----------------------------------------------------------------------------
proc rhsv2(uc, uv, ref Rv)
{
   const a = 1.0/3.0, b = -2.0, c = 1.0, d = 2.0/3.0;

   forall i in D
   {
      const J = jacobian(uv[i]);
      // cell center values
      const uim1 = (6.0 * uc[i-1] - uv[i-1] - uv[i]) / 4.0;
      const ui   = (6.0 * uc[i]   - uv[i+1] - uv[i]) / 4.0;
      var ux : real;
      if J > 0.0
      {
         // backward difference
         ux = (a * uv[i-1] + b * uim1 + c * uv[i] + d * ui) / dx;
      }
      else
      {
         // forward difference
         ux = -(d * uim1 + c * uv[i] + b * ui + a * uv[i+1]) / dx;
      }
      Rv[i] = J * ux;
   }
}

//-----------------------------------------------------------------------------
// RHS at cell centers
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
         if diff == 1 then
            rhsv1(uc, uv, Rv);
         else
            rhsv2(uc, uv, Rv);

         rhsc(uv, Rc);

         forall i in D
         {
            uc[i] = ark[rk] * uc0[i] + brk[rk] * (uc[i] - dt * Rc[i]);
            uv[i] = ark[rk] * uv0[i] + brk[rk] * (uv[i] - dt * Rv[i]);
         }
         uc.updateFluff();
         uv.updateFluff();
      }
      t = t + dt;
      writeln("t, dt = ", t, " ", dt);
   }

   output(xc, uc, "uc.txt");
   output(xv, uv, "uv.txt");
}
