/*
 Compile as
     g++ main.cc
 Run the program by specifying input file
     ./a.out test.in
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cassert>

#define SIGN(a) (((a)<0) ? -1:1)

const double arks[] = {0.0, 3.0/4.0, 1.0/3.0};
const double brks[] = {1.0, 1.0/4.0, 2.0/3.0};

enum FluxScheme { roe, godunov };
enum ReconstructScheme { first, muscl_minmod, muscl_vanleer };

using namespace std;

//------------------------------------------------------------------------------
// Minmod of three numbers
//------------------------------------------------------------------------------
double minmod (const double& a,
               const double& b,
               const double& c)
{
   double result;
   
   if (a*b > 0.0 && b*c > 0.0)
   {
      result  = min( min(fabs(a), fabs(b)), fabs(c) );
      result *= SIGN(a);
   }
   else
      result = 0.0;
   
   return result;
   
}

//------------------------------------------------------------------------------
// vanleer limiter
//------------------------------------------------------------------------------
double vanleer (const double& a,
                const double& b)
{
   double du;

   if(fabs(a*b) > 0.0)
      du = (SIGN(a)+SIGN(b))*fabs(a)*fabs(b)/(fabs(a) + fabs(b));
   else
      du = 0.0;

   return du;
}
//------------------------------------------------------------------------------
// Reconstruct left state of right face
//------------------------------------------------------------------------------
double muscl (const double& ul,
              const double& uc,
              const double& ur)
{
   static const double beta = 2.0;
   
   double dul = uc - ul;
   double dur = ur - uc;
   double duc = (ur - ul)/2.0;
   double result = uc + 0.5 * minmod (beta*dul, duc, beta*dur);
   
   return result;
}

//------------------------------------------------------------------------------
// Flux of burgers equation
//------------------------------------------------------------------------------
double phys_flux(const double u)
{
   return 0.5*u*u;
}

//------------------------------------------------------------------------------
// Main class of the problem
//------------------------------------------------------------------------------
class FVProblem
{
   public:
      FVProblem (char param_file[]);
      void run ();
   
   private:
      void read_parameters (char param_file[]);
      void make_grid_and_dofs ();
      void initial_condition ();
      void compute_dt ();
      void reconstruct (const unsigned int face,
                        double& left,
                        double& right) const;
      void roe_flux (const double,
                     const double,
                           double&) const;
      void num_flux (const double,
                     const double,
                     double&) const;
      void compute_residual ();
      void update_solution (const unsigned int rk);
      void output ();
   
      //Parameters
      double u_left, u_right;
      unsigned int n_cell;
      unsigned int n_face;
      double xmin;
      double xmax;
      double xmid;
      double dx;
      double dt;
      double cfl;
      double final_time;
      FluxScheme flux_scheme;
   
      // Grid
      vector<double> xc;
      vector<double> xf;
   
      vector<double> residual;
      vector<double> conserved;
      vector<double> conserved_old;
};

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
FVProblem::FVProblem (char param_file[])
{
   read_parameters (param_file);
   
}
//------------------------------------------------------------------------------
// Read parameters from file
//------------------------------------------------------------------------------
void FVProblem::read_parameters (char param_file[])
{
   cout << "Reading parameters from file " << param_file << endl;
   
   ifstream fp(param_file);
   assert (fp.is_open());
   
   string param, input;
   
   fp >> param >> n_cell;
   fp >> param >> xmin;
   fp >> param >> xmax;
   fp >> param >> xmid;
   fp >> param >> u_left;
   fp >> param >> u_right;
   fp >> param >> cfl;
   fp >> param >> final_time;
   fp >> param >> input;
   if(input == "roe")
      flux_scheme = roe;
   else
   {
      cout << "Unknown flux scheme\n";
      exit(0);
   }
   
   fp.close();
}

//------------------------------------------------------------------------------
// Allocate memory for grid and create grid
//------------------------------------------------------------------------------
void FVProblem::make_grid_and_dofs ()
{
   cout << "Making grid and allocating memory ...\n";
   
   n_face = n_cell + 1;
   dx = (xmax - xmin) / n_cell;
   xc.resize (n_cell);
   xf.resize (n_face);
   
   // Make grid
   for(unsigned int i=0; i<n_face; ++i)
      xf[i] = xmin + i * dx;
   for(unsigned int i=0; i<n_cell; ++i)
      xc[i] = 0.5 * (xf[i] + xf[i+1]);
   
   residual.resize (n_cell);
   conserved.resize (n_cell);
   conserved_old.resize (n_cell);
}

//------------------------------------------------------------------------------
// Set initial condition
//------------------------------------------------------------------------------
void FVProblem::initial_condition ()
{
   cout << "Setting initial conditions ...\n";
   
   // Set initial condition
   for(unsigned int i=0; i<n_cell; ++i)
   {
      if(xf[i+1] <= xmid)
         conserved[i] = u_left;
      else
         conserved[i] = u_right;
   }
}

//------------------------------------------------------------------------------
// Compute time step
//------------------------------------------------------------------------------
void FVProblem::compute_dt ()
{
   dt = 1.0e20;
   for(unsigned int i=0; i<n_cell; ++i)
   {
      double speed = fabs(conserved[i]) + 1.0e-14;
      dt = min (dt, dx/speed);
   }
   dt *= cfl;
}

//------------------------------------------------------------------------------
// Reconstruct left/right state at a face
//------------------------------------------------------------------------------
void FVProblem::reconstruct (const unsigned int face,
                             double& left,
                             double& right) const
{
   if(face==1)
   {
      left = conserved[0];
      right = muscl (conserved[2], conserved[1], conserved[0]);
   }
   else if(face==n_face-2)
   {
      left  = muscl (conserved[n_cell-3], conserved[n_cell-2],
                     conserved[n_cell-1]);
      right = conserved[n_cell-1];
   }
   else
   {
      left  = muscl (conserved[face-2], conserved[face-1], conserved[face]);
      right = muscl (conserved[face+1], conserved[face], conserved[face-1]);
   }
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void FVProblem::roe_flux(const double left,
                         const double right,
                               double& flux) const
{
   double a = fabs( 0.5*(left + right) );
   flux = 0.5*( phys_flux(left) + phys_flux(right) ) - 0.5 * a * (right - left);
}

//------------------------------------------------------------------------------
// Numerical flux function
//------------------------------------------------------------------------------
void FVProblem::num_flux(const double left,
                         const double right,
                         double&       flux) const
{
   switch(flux_scheme)
   {
      case roe:
         roe_flux (left, right, flux);
         break;
         
      default:
         cout << "Unknown flux scheme specified\n";
         exit(0);
   }
}

//------------------------------------------------------------------------------
// Compute finite volume residual
//------------------------------------------------------------------------------
void FVProblem::compute_residual ()
{
   for(unsigned int i=0; i<n_cell; ++i)
      residual[i] = 0.0;
   
   double flux, left, right;
   
   // Flux through left boundary
   num_flux (u_left, conserved[0], flux);
   residual[0] -= flux;

   // Flux through interior faces
   for(unsigned int i=1; i<n_face-1; ++i)
   {
      reconstruct (i, left, right);
      num_flux (left, right, flux);
      
      residual[i-1] += flux;
      residual[i]   -= flux;
   }

   // Flux through right boundary
   num_flux (conserved[n_cell-1], u_right, flux);
   residual[n_cell-1] += flux;
}

//------------------------------------------------------------------------------
// Compute finite volume residual
//------------------------------------------------------------------------------
void FVProblem::update_solution (const unsigned int rk)
{
   for(unsigned int i=0; i<n_cell; ++i)
         conserved[i] = arks[rk] * conserved_old[i] +
            brks[rk] * (conserved[i] - (dt/dx) * residual[i]);
   
}

//------------------------------------------------------------------------------
// Save solution to file
//------------------------------------------------------------------------------
void FVProblem::output ()
{
   cout << "Saving solution to sol.dat\n";
   
   ofstream fo("sol.dat");
   for(unsigned int i=0; i<n_cell; ++i)
      fo << xc[i] << " "
         << conserved[i] << endl;
   fo.close ();
}

//------------------------------------------------------------------------------
// Start the computations
//------------------------------------------------------------------------------
void FVProblem::run ()
{
   make_grid_and_dofs ();
   initial_condition ();

   double time = 0.0;
   unsigned int iter = 0;
   while (time < final_time)
   {
      conserved_old = conserved;
      compute_dt ();
      if(time+dt > final_time) dt = final_time - time;
      for(unsigned int rk=0; rk<3; ++rk)
      {
         compute_residual ();
         update_solution (rk);
      }
      time += dt;
      ++iter;
      if(iter % 1000 == 0) 
      cout << "Iter = " << iter << " Time = " << time << endl;
      if(iter % 1000 == 0) output ();
   }
   
   output ();
}

//------------------------------------------------------------------------------
int main (int argv, char* argc[])
{
   assert(argv == 2);
   FVProblem fv_problem (argc[1]);
   fv_problem.run ();
   
   return 0;
}
