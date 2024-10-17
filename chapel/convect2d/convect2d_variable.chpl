/* Solves u_t + div(vel*u) = 0 using finite volume scheme, 
   where vel = (-y,x)
   weno5 reconstruction and periodic boundary conditions. 
   You need to use StencilDist which will be available in v1.14
   of Chapel. You can specify some command line options, e.g.,
      ./convect2d_variable --nx=100 --ny=100 --Tf=5.0 --cfl=0.8 --si=100
   See below for explanation of nx, ny, Tf, cfl, si. 
   Solution is saved in Tecplot format which can also be opened
   in VisIt.
      visit -o sol*.tec
*/
use IO;
use StencilDist;
use Math;

config const nx = 100,   // number of cells in x direction
             ny = 100,   // number of cells in y direction
             Tf = 10.0, // Time of simulation
             cfl= 0.9,  // cfl number
             si = 100;  // iteration interval to save solution
const xmin = -1.0, xmax = 1.0,
      ymin = -1.0, ymax = 1.0;
const ark : [1..3] real = (0.0, 3.0/4.0, 1.0/3.0);
const brk : [1..3] real = (1.0, 1.0/4.0, 2.0/3.0);
var dx, dy : real;

// (x,y) : coordinates
// vel   : velocity at (x,y)
proc advection_velocity(x:real, y:real, ref vel:[1..2] real)
{
   vel[1] = -y;
   vel[2] =  x;
}

// fv weno5 reconstruction, gives left state at interface
// between u0 and up1
proc weno5(um2:real, um1:real, u0:real, up1:real, up2:real) : real
{
   const eps = 1.0e-6;
   const gamma1=1.0/10.0, gamma2=3.0/5.0, gamma3=3.0/10.0;

   const beta1 = (13.0/12.0)*(um2 - 2.0*um1 + u0)**2 +
                 (1.0/4.0)*(um2 - 4.0*um1 + 3.0*u0)**2;
   const beta2 = (13.0/12.0)*(um1 - 2.0*u0 + up1)**2 +
                 (1.0/4.0)*(um1 - up1)**2;
   const beta3 = (13.0/12.0)*(u0 - 2.0*up1 + up2)**2 +
                 (1.0/4.0)*(3.0*u0 - 4.0*up1 + up2)**2;

   const w1 = gamma1 / (eps+beta1)**2;
   const w2 = gamma2 / (eps+beta2)**2;
   const w3 = gamma3 / (eps+beta3)**2;

   const u1 = (1.0/3.0)*um2 - (7.0/6.0)*um1 + (11.0/6.0)*u0;
   const u2 = -(1.0/6.0)*um1 + (5.0/6.0)*u0 + (1.0/3.0)*up1;
   const u3 = (1.0/3.0)*u0 + (5.0/6.0)*up1 - (1.0/6.0)*up2;

   return (w1 * u1 + w2 * u2 + w3 * u3)/(w1 + w2 + w3);
}

// upwind flux
// (lx,ly): unit normal vector to face
// vel    : velocity vector
// ul,ur  : left, right states
// flux   : numerical flux
proc numerical_flux(lx:real, ly:real, vel:[1..2] real, 
                    ul:real, ur:real, ref flux:real)
{
   const vn = vel[1]*lx + vel[2]*ly; // normal velocity
   if(vn > 0.0)
   {
      flux = vn * ul;
   }
   else
   {
      flux = vn * ur;
   }
}

// Save solution to file
proc savesol(t : real, u : [?D] real, c : int) : int
{
  // construct filename with counter c
  if c > 999 then halt("Filename counter too large !!!");
  const filename = "sol%03i.tec".format(c);

  var fw = open(filename, ioMode.cw).writer(locking=false);
  fw.writeln("TITLE = \"u_t + u_x + u_y = 0\"");
  fw.writeln("VARIABLES = x, y, sol");
  fw.writeln("ZONE STRANDID=1, SOLUTIONTIME=",t,", I=",nx,", J=",ny,", DATAPACKING=POINT");
  for (j,i) in D
  {
    const x = xmin + (i-1)*dx + 0.5*dx,
          y = ymin + (j-1)*dy + 0.5*dy;
    fw.writeln(x,"  ",y,"  ",u[i,j]);
  }
  fw.close();
  return c+1;
}

// main function
proc main()
{
  dx = (xmax-xmin)/nx;
  dy = (ymax-ymin)/ny;

  writeln("Grid size is ",nx," x ",ny);
  writeln("dx, dy =", dx, dy);

  const D  = {1..nx, 1..ny},     // nx * ny     cells
        Dx = {1..(nx+1),1..ny},  // (nx+1) * ny vertical faces
        Dy = {1..nx,1..(ny+1)};  // nx * (ny+1) horizontal faces
  param rank = D.rank;

  // problem space
  var halo: rank*int = (3,3);
  const PSpace = D dmapped new stencilDist(D, fluff=halo, periodic=true);

  var u : [PSpace] real;

  // Set initial condition
  forall (i,j) in D
  {
    const x = xmin + (i-1)*dx + 0.5*dx,
          y = ymin + (j-1)*dy + 0.5*dy;
    u[i,j] = 1.0 + exp(-100.0*((x-0.5)**2 + y**2));
  }
  u.updateFluff();

  var c  = 0;     // counter to save solution
  c = savesol(0.0, u, c);

  // Compute dt
  var dt = 1.0e20;
  for (i,j) in D
  {
    const x = xmin + (i-1)*dx + 0.5*dx,
          y = ymin + (j-1)*dy + 0.5*dy; // cell center
    var vel : [1..2] real;
    advection_velocity(x, y, vel);
    dt = min(dt, 1.0/(abs(vel[1])/dx + abs(vel[2])/dy) + 1.0e-20);
  }
  dt *= cfl;

  var u0 : [PSpace] real; // old solution
  var res: [PSpace] real; // residual

  var t  = 0.0;   // time counter
  var it = 0;     // iteration counter
  var lam = dt/(dx*dy);
  while t < Tf
  {
    // Adjust dt so we exactly reach Tf
    if t+dt > Tf 
    {
      dt  = Tf - t;
      lam = dt/(dx*dy);
    }

    u0 = u;

    // 3-stage RK scheme
    // du/dt + res(u) = 0
    for rk in 1..3
    {
      res = 0.0;
      // x fluxes
      forall (i,j) in Dx
      {
        // velocity at face mid-point
        const xf = xmin + (i-1)*dx,
              yf = ymin + (j-1)*dy + 0.5*dy;
        var vel : [1..2] real;
        advection_velocity(xf, yf, vel);
        // reconstruct left/right states at face between (i-1,j) and (i,j)
        const ul = weno5(u[i-3,j],u[i-2,j],u[i-1,j],u[i,j],u[i+1,j]);
        const ur = weno5(u[i+2,j],u[i+1,j],u[i,j],u[i-1,j],u[i-2,j]);
        // compute numerical flux
        var flux : real;
        numerical_flux(1.0, 0.0, vel, ul, ur, flux);
        res[i-1,j] += flux * dy;
        res[i,j]   -= flux * dy;
      }
      // y fluxes
      forall (i,j) in Dy
      {
        // velocity at face mid-point
        const xf = xmin + (i-1)*dx + 0.5*dx,
              yf = ymin + (j-1)*dy;
        var vel : [1..2] real;
        advection_velocity(xf, yf, vel);
        // reconstruct left/right states at face between (i,j-1) and (i,j)
        const ul = weno5(u[i,j-3],u[i,j-2],u[i,j-1],u[i,j],u[i,j+1]);
        const ur = weno5(u[i,j+2],u[i,j+1],u[i,j],u[i,j-1],u[i,j-2]);
        // compute numerical flux
        var flux : real;
        numerical_flux(0.0, 1.0, vel, ul, ur, flux);
        res[i,j-1] += flux * dx;
        res[i,j]   -= flux * dx;
      }
      u = ark[rk] * u0 + brk[rk] * (u - lam * res);
      u.updateFluff();
    }

    t  += dt;
    it += 1;
    writeln(it,"  ",t);
    if it%si == 0 then c = savesol(t, u, c);
  }
  // Save solution at final time
  c = savesol(t, u, c);
}
