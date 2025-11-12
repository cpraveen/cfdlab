use Math;
use Poisson;
use VTK;

config const nx     = 128,
             ny     = 128,
             levels = 7,
             rtol   = 1.0e-6,
             niter  = 1000,
             nsmooth = 2,
             method  = "mg";

const Dx = 0..nx;
const Dy = 0..ny;
const D = {Dx, Dy};
const xmin = 0.0, xmax = 1.0;
const ymin = 0.0, ymax = 1.0;
const dx = (xmax - xmin) / nx;
const dy = (ymax - ymin) / ny;

proc rhs(x,y) do return 2 * (2*pi)**2 * sin(2*pi*x) * sin(2*pi*y);

proc main()
{
   // Make 1d grids
   const x = [i in Dx] xmin + i*dx;
   const y = [j in Dy] ymin + j*dy;

   var v, f : [D] real;

   // Set rhs, only inner points needed
   forall (i,j) in D.expand(-1) do
      f[i,j] = rhs(x[i], y[j]);

   if method == "mg" then
      multigrid(v, f, dx, dy, levels, rtol, niter, nsmooth);
   else
      sor(v, f, dx, dy, rtol, niter);

   write_vtk(x,y,"sol",v,"sol.vtk");
}
