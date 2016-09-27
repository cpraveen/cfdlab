/* Solves u_t + u_x + u_y = 0 using finite volume scheme, 
   weno5 reconstruction and periodic boundary conditions. 
   You need to use StencilDist which will be available in v1.14
   of Chapel. You can specify some command line options, e.g.,
      ./a.out --n=100 --Tf=5.0 --cfl=0.4 --si=100
   See below for explanation of n, Tf, cfl, si. 
   Solution is saved in Tecplot format which can also be opened
   in VisIt.
*/
use StencilDist;

config const n  = 50;   // number of cells in each direction
config const Tf = 10.0; // Time of simulation
config const cfl= 0.4;  // cfl number
config const si = 100;  // iteration interval to save solution
const xmin = 0.0, xmax = 1.0;
const ymin = 0.0, ymax = 1.0;
const ark : [1..3] real = [0.0, 3.0/4.0, 1.0/3.0];
const brk = 1.0 - ark;
var dx, dy : real;

// fv weno5 reconstruction, gives left state at interface
// b/w u0 and u1
proc weno5(um2:real, um1:real, u0:real, up1:real, up2:real): real
{
   const eps = 1.0e-6;
   const gamma1=1.0/10.0, gamma2=3.0/5.0, gamma3=3.0/10.0;
   var beta1, beta2, beta3 : real;
   var u1, u2, u3 : real;
   var w1, w2, w3 : real;

   beta1 = (13.0/12.0)*(um2 - 2.0*um1 + u0)**2 +
           (1.0/4.0)*(um2 - 4.0*um1 + 3.0*u0)**2;
   beta2 = (13.0/12.0)*(um1 - 2.0*u0 + up1)**2 +
           (1.0/4.0)*(um1 - up1)**2;
   beta3 = (13.0/12.0)*(u0 - 2.0*up1 + up2)**2 +
           (1.0/4.0)*(3.0*u0 - 4.0*up1 + up2)**2;

   w1 = gamma1 / (eps+beta1)**2;
   w2 = gamma2 / (eps+beta2)**2;
   w3 = gamma3 / (eps+beta3)**2;

   u1 = (1.0/3.0)*um2 - (7.0/6.0)*um1 + (11.0/6.0)*u0;
   u2 = -(1.0/6.0)*um1 + (5.0/6.0)*u0 + (1.0/3.0)*up1;
   u3 = (1.0/3.0)*u0 + (5.0/6.0)*up1 - (1.0/6.0)*up2;

   return (w1 * u1 + w2 * u2 + w3 * u3)/(w1 + w2 + w3);
}

// Save solution to file
proc savesol(t : real, u : [?D] real, c : int) : int
{
  // construct filename with counter c
  var filename : string;
  if c < 9
  {
    filename = "sol000"+c:string+".tec";
  }
  else if c < 99 
  {
    filename = "sol00"+c:string+".tec";
  }
  else if c < 999 
  {
    filename = "sol0"+c:string+".tec";
  }
  else
  {
    writeln("Filename counter too large !!!");
    exit();
  }

  var fw = open(filename, iomode.cw).writer();
  fw.writeln("TITLE = \"u_t + u_x + u_y = 0\"");
  fw.writeln("VARIABLES = x, y, sol");
  fw.writeln("ZONE STRANDID=1, SOLUTIONTIME=",t,", I=",n,", J=",n,", DATAPACKING=POINT");
  for j in 1..n
  {
    for i in 1..n
    {
      var x = xmin + (i-1)*dx + 0.5*dx;
      var y = ymin + (j-1)*dy + 0.5*dy;
      fw.writeln(x,"  ",y,"  ",u[i,j]);
    }
  }
  return c+1;
}

// main function
proc main()
{
  dx = (xmax-xmin)/n;
  dy = (ymax-ymin)/n;

  writeln("Grid size is ",n," x ",n);
  writeln("dx, dy =", dx, dy);

  const dom = {1..n, 1..n};
  param rank = dom.rank;

  var halo1, halo3 : rank*int;
  for i in 1..rank do halo1(i) = 1;
  for i in 1..rank do halo3(i) = 3;
  const  Space = dom dmapped Stencil(dom, fluff=halo1, periodic=false);
  const PSpace = dom dmapped Stencil(dom, fluff=halo3, periodic=true);

  var u : [PSpace] real;

  // Set initial condition
  forall (i,j) in PSpace
  {
    var x = xmin + (i-1)*dx + 0.5*dx;
    var y = ymin + (j-1)*dy + 0.5*dy;
    u[i,j] = sin(2*pi*x) * sin(2*pi*y);
  }
  u.updateFluff();

  var c  = 0;     // counter to save solution
  c = savesol(0.0, u, c);

  var u0 : [PSpace] real;
  var res: [PSpace] real;

  const umax = sqrt(2.0); // wave speed used in dt computation
  var dt = cfl * dx / umax;
  var t  = 0.0;   // time counter
  var it = 0;     // iteration counter
  const lam = dt/(dx*dy);
  while t < Tf
  {
    // Adjust dt so we exactly reach Tf
    if t+dt > Tf then dt = Tf - t;

    u0 = u;

    // 3-stage RK scheme
    for rk in 1..3
    {
      res = 0.0;
      // x fluxes
      forall (i,j) in {1..(n+1),1..n}
      {
        var ul = weno5(u[i-3,j],u[i-2,j],u[i-1,j],u[i,j],u[i+1,j]);
        res[i-1,j] += ul * dy;
        res[i,j]   -= ul * dy;
      }
      // y fluxes
      forall (i,j) in {1..n,1..(n+1)}
      {
        var ul = weno5(u[i,j-3],u[i,j-2],u[i,j-1],u[i,j],u[i,j+1]);
        res[i,j-1] += ul * dx;
        res[i,j]   -= ul * dx;
      }
      u = ark[rk] * u0 + brk[rk] * (u - lam * res);
      u.updateFluff();
    }

    t  += dt;
    it += 1;
    writeln(it,"  ",t);

    if it%si == 0 then c = savesol(t, u, c);
  }

  c = savesol(t, u, c);
}
