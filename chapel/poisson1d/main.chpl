// Solve -u'' = f

use Math;
use IO;
use Poisson;

config const N       = 128,     // grid of N+1 points
             levels  = 6,       // number of levels in V-cycle
             rtol    = 1.0e-6,  // relative tolerance for vcycle
             niter   = 5000,    // max number of v-cycles
             nsmooth = 2,       // number of jacobi iterations
             method  = "mg";    // mg, sor

const xmin = 0.0;
const xmax = 1.0;

proc exact(x) do return x + sin(2*pi*x);
proc rhs(x) do return (2*pi)**2 * sin(2*pi*x);

//------------------------------------------------------------------------------
proc main()
{
   const h = (xmax - xmin) / N;
   const D = {1..N+1};

   var x, v, f, vnew, r : [D] real;
   forall i in D
   {
      x[i] = xmin + (i-1)*h;
      f[i] = rhs(x[i]);
   }

   // Initial v must have correct bc filled in
   v[1]   = exact(xmin);
   v[N+1] = exact(xmax);

   if method == "mg"
   {
      assert(N % 2**levels == 0);
      multigrid(v, f, h, rtol, niter, nsmooth, levels);
   }
   else
   {
      sor(v, f, h, rtol, niter);
   }

   const filename = "sol.txt";
   var fw = open(filename, ioMode.cw).writer(locking=false);
   for i in D
   {
      fw.writef("%14.8er %14.8er\n", x[i], v[i]);
   }
   fw.close();
}
