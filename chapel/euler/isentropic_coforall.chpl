/* Solves isentropic vortex for 2D Euler equation using finite volume scheme with the
   Lax Friedrich flux, WENO5 reconstruction and periodic boundary conditions.
   You need to use StencilDist which will be available in v1.14
   of Chapel. You can specify some command line options, e.g.,
      ./isentropic --nx=100 --ny=100 --Tf=20.0 --cfl=0.4 --si=100
   See below for explanation of n, Tf, cfl, si.
   Solution is saved in Tecplot format which can also be opened
   in VisIt.
      visit -o sol*.tec
   Original Author: Deep Ray, 8 Oct 2016
*/
use IO;
use StencilDist;
use Math;

config const nx    = 50,           // number of cells in x direction
             ny    = 50,           // number of cells in y direction
             cfl   = 0.45,         // cfl number
             si    = 100,          // iteration interval to save solution
             M     = 0.5,          // Farfield Mach number of flow
             alpha = 45.0,         // Angle of advection (in degrees)
             Tf    = sqrt(2)*10/M; // Time of simulation
const xmin = -5.0, xmax = 5.0,
      ymin = -5.0, ymax = 5.0;
const ark : [1..3] real = (0.0, 3.0/4.0, 1.0/3.0);
const brk : [1..3] real = (1.0, 1.0/4.0, 2.0/3.0);
const gamma = 1.4;
const gas_const = 1.0;
var dx, dy : real;
const nvar = 4;
type cType = [1..nvar] real;

// fv weno5 reconstruction, gives left state at interface
// between u0 and up1
proc weno5(Um2:cType, Um1:cType, U0:cType, Up1:cType, Up2:cType): cType
{
   const eps = 1.0e-6;
   const gamma1=1.0/10.0, gamma2=3.0/5.0, gamma3=3.0/10.0;
   var U : cType;

   for i in 1..nvar
   {
      const um2 = Um2[i];
      const um1 = Um1[i];
      const u0  = U0[i];
      const up1 = Up1[i];
      const up2 = Up2[i];

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

      U[i] =  (w1 * u1 + w2 * u2 + w3 * u3)/(w1 + w2 + w3);
   }

   return U;
}

// Save solution to file
proc savesol(t : real, U : [?D] cType, c : int) : int
{
  // construct filename with counter c
  if c > 999 then halt("Filename counter too large !!!");
  const filename = "sol%03i.tec".format(c);

  var fw = open(filename, ioMode.cw).writer();
  fw.writeln("TITLE = \"u_t + u_x + u_y = 0\"");
  fw.writeln("VARIABLES = x, y, rho, u, v, p");
  fw.writeln("ZONE STRANDID=1, SOLUTIONTIME=",t,", I=",nx,", J=",ny,", DATAPACKING=POINT");
  for (j,i) in D
  {
    const x = xmin + (i-1)*dx + 0.5*dx,
          y = ymin + (j-1)*dy + 0.5*dy;
    const Prim = con2prim(U[i,j]);
    fw.writeln(x,"  ",y,"  ",Prim[1],"  ",Prim[2],"  ",Prim[3],"  ",Prim[4]);
  }
  fw.close();
  return c+1;
}

// Conserved to primitive variables
proc con2prim(Con : cType) : cType
{
  var Prim : cType;
  Prim[1] = Con[1];
  Prim[2] = Con[2]/Con[1];
  Prim[3] = Con[3]/Con[1];
  Prim[4] = (Con[4] - 0.5*Prim[1]*(Prim[2]**2 + Prim[3]**2))*(gamma-1.0);
  return Prim;
}

// Primitive to conserved variables
proc prim2con(Prim : cType) :  cType
{
  var Con : cType;
  Con[1] = Prim[1];
  Con[2] = Prim[1]*Prim[2];
  Con[3] = Prim[1]*Prim[3];
  Con[4] = 0.5*Prim[1]*(Prim[2]**2 + Prim[3]**2) + Prim[4]/(gamma-1.0);
  return Con;
}

// Compute maximum eigenvalue in direction (n1,n2)
proc maxeigval(Con : cType, n1 : real, n2 : real) : real
{
   const Prim  = con2prim(Con);
   const u     = abs(Prim[2]*n1 + Prim[3]*n2);
   const a     = sqrt(gamma*Prim[4]/Prim[1]);
   return u + a;
}

// Compute local timestep
proc dt_local(Con : cType) : real
{
   const Prim  = con2prim(Con);
   const u     = sqrt(Prim[2]**2 + Prim[3]**2);
   const a     = sqrt(gamma*Prim[4]/Prim[1]);
   const eigen = u + a;
   return min(dx,dy)/eigen;
}

// Rusanov flux
proc Flux(Ul : cType, Ur : cType, n1 : real, n2 : real) : cType
{
  const U = 0.5*(Ul + Ur);
  const lam = maxeigval(U, n1, n2 );
  return avg_flux(Ul,Ur,n1,n2) - 0.5*lam*(Ur-Ul);
}

// Simple average flux
proc avg_flux(Ul : cType, Ur : cType, n1 : real, n2 : real) : cType
{
  const Pl = con2prim(Ul);
  const Pr = con2prim(Ur);

  var fluxl, fluxr: cType;

  fluxl[1] = Ul[2]*n1 + Ul[3]*n2;
  fluxl[2] = Pl[4]*n1 + Pl[2]*fluxl[1];
  fluxl[3] = Pl[4]*n2 + Pl[3]*fluxl[1];
  fluxl[4] = (Ul[4]+Pl[4])*(Pl[2]*n1 + Pl[3]*n2);

  fluxr[1] = Ur[2]*n1 + Ur[3]*n2;
  fluxr[2] = Pr[4]*n1 + Pr[2]*fluxr[1];
  fluxr[3] = Pr[4]*n2 + Pr[3]*fluxr[1];
  fluxr[4] = (Ur[4]+Pr[4])*(Pr[2]*n1 + Pr[3]*n2);

  return 0.5*(fluxl+fluxr);
}

// main function
proc main()
{
  dx = (xmax-xmin)/nx;
  dy = (ymax-ymin)/ny;

  writeln("Grid size is ",nx," x ",ny);
  writeln("dx, dy = ", dx, "  ", dy);

  const D  = {1..nx  , 1..ny},
        Dx = {1..nx+1, 1..ny},
        Dy = {1..nx  , 1..ny+1},
        Dr = {0..nx+1, 0..ny+1};
  param rank = D.rank;

  var halo: rank*int = (3,3);
  const PSpace = stencilDist.createDomain(D, fluff=halo, periodic=true);

  var U, U0 : [PSpace] cType;
  var res   : [Dr] cType;

  // Set initial condition
  const beta = 5.0;
  forall (i,j) in D
  {
    const x = xmin + (i-1)*dx + 0.5*dx,
          y = ymin + (j-1)*dy + 0.5*dy;
    const r2 = x**2 + y**2;
    var Prim : cType;
    const T = 1.0 - (gamma-1.0)*(beta**2)/(8.0*gamma*pi*pi)*exp(1-r2);
    Prim[1] = T**(1.0/(gamma-1.0));
    Prim[2] = M*cos(alpha*pi/180.0) - beta/(2.0*pi)*y*exp(0.5*(1.0-r2));
    Prim[3] = M*sin(alpha*pi/180.0) + beta/(2.0*pi)*x*exp(0.5*(1.0-r2));
    Prim[4] = Prim[1]**gamma;
    U[i,j] = prim2con(Prim);
  }
  U.updateFluff();

  var c  = 0;     // counter to save solution
  c = savesol(0.0, U, c);
  var dt, lam : real;
  var t  = 0.0;   // time counter
  var it = 0;     // iteration counter
  while t < Tf
  {
    dt = min reduce dt_local(U);
    dt *= cfl;
    // Adjust dt so we exactly reach Tf
    if t+dt > Tf then dt = Tf - t;
    lam = dt/(dx*dy);

    U0 = U;

    // 3-stage RK scheme
    for rk in 1..3
    {
      forall (i,j) in Dr
      {
        res[i,j] = 0.0;
      }

      // x fluxes
      coforall j in 1..ny
      {
        for i in 1..nx+1
        {
          const Ul = weno5(U[i-3,j],U[i-2,j],U[i-1,j],U[i,j],U[i+1,j]);
          const Ur = weno5(U[i+2,j],U[i+1,j],U[i,j],U[i-1,j],U[i-2,j]);
          const flux = dy * Flux(Ul,Ur,1.0,0.0);
          res[i-1,j] += flux;
          res[i  ,j] -= flux;
        }
      }

      // y fluxes
      coforall i in 1..nx
      {
        for j in 1..ny+1
        {
          const Ul = weno5(U[i,j-3],U[i,j-2],U[i,j-1],U[i,j],U[i,j+1]);
          const Ur = weno5(U[i,j+2],U[i,j+1],U[i,j],U[i,j-1],U[i,j-2]);
          const flux = dx * Flux(Ul,Ur,0.0,1.0);
          res[i,j-1] += flux;
          res[i,j  ] -= flux;
        }
      }

      forall (i,j) in D
      {
         U[i,j]  = ark[rk]*U0[i,j] + brk[rk]*(U[i,j] - lam*res[i,j]);
      }

      U.updateFluff();
    }

    t  += dt;
    it += 1;
    writeln(it,"  ",t);
    if it%si == 0 then c = savesol(t, U, c);
  }
  // Save solution at final time
  c = savesol(t, U, c);
}
