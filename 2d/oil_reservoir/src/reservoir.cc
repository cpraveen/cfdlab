#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cassert>
#include "material.h"
#include "pressure.h"
#include "reservoir.h"

#define SIGN(a) (((a)<0) ? -1:1)

extern int n_interior_min;
const double beta = 2.0; // factor in minmod limiter

using namespace std;

//------------------------------------------------------------------------------
// Minmod limiter
// minmod ( 2(u0-ul), (ur-ul)/2, 2(ur-u0) )
//------------------------------------------------------------------------------
double minmod (const double& ul, const double& u0, const double& ur)
{
   double result;

   double db = u0 - ul;         // backward difference
   double df = ur - u0;         // forward difference
   double dc = 0.5 * (ur - ul); // central difference

   if (db*dc > 0.0 && dc*df > 0.0)
   {
      result = min( min(fabs(beta*db), fabs(dc)), fabs(beta*df) );
      result *= SIGN(db);
   }
   else
      result = 0.0;

   return result;
}

//------------------------------------------------------------------------------
// Find left state at interface between (il,jl) and (ir,jr)
// (ill,jll) is to the left of (il,jl)
//------------------------------------------------------------------------------
vector<double> ReservoirProblem::reconstruct
       (
       const unsigned int& ill,
       const unsigned int& jll,
       const unsigned int& il,
       const unsigned int& jl,
       const unsigned int& ir,
       const unsigned int& jr
       ) const
{
   if (order==1)
      return reconstruct1(ill, jll, il, jl, ir, jr);
   else
      return reconstruct2(ill, jll, il, jl, ir, jr);
}

//------------------------------------------------------------------------------
// First order reconstruction
//------------------------------------------------------------------------------
vector<double> ReservoirProblem::reconstruct1
       (
       const unsigned int& ill,
       const unsigned int& jll,
       const unsigned int& il,
       const unsigned int& jl,
       const unsigned int& ir,
       const unsigned int& jr
       ) const
{
   vector<double> state(3);

   state[0] = saturation    (il, jl);
   state[1] = concentration (il, jl); 
   state[2] = permeability  (il, jl); 

   return state;
}

//------------------------------------------------------------------------------
// Second order reconstruction
//------------------------------------------------------------------------------
vector<double> ReservoirProblem::reconstruct2
       (
       const unsigned int& ill,
       const unsigned int& jll,
       const unsigned int& il,
       const unsigned int& jl,
       const unsigned int& ir,
       const unsigned int& jr
       ) const
{
   vector<double> state(3);
   double ds, dp;

   // saturation
   ds = minmod (saturation(ill,jll), 
                saturation(il,jl), 
                saturation(ir,jr));
   state[0] = saturation (il,jl) + 0.5 * ds;

   // concentration
   // We reconstruct b = s*c
   double bll = saturation(ill,jll) * concentration(ill,jll);
   double bl  = saturation(il ,jl ) * concentration(il ,jl );
   double br  = saturation(ir ,jr ) * concentration(ir ,jr );
   double db  = minmod (bll, bl, br);
   state[1]   = (bl + 0.5 * db);
   if (state[0] > SZERO)
      state[1] /= state[0];
   else
      state[1] = 0.0;

   // permeability
   dp = minmod (permeability(ill,jll), 
                permeability(il,jl), 
                permeability(ir,jr));
   state[2] = permeability (il, jl) + 0.5 * dp;

   return state;
}

//------------------------------------------------------------------------------
// computes total darcy velocity at interface b/w
// (ileft,jleft) and (iright,jright)
//------------------------------------------------------------------------------
double ReservoirProblem::darcy_velocity
       (
       const unsigned int& ileft,
       const unsigned int& jleft,
       const unsigned int& iright,
       const unsigned int& jright,
       const double& g
       )
{
   double s_left  = saturation    (ileft,  jleft);
   double c_left  = concentration (ileft,  jleft);
   double p_left  = pressure      (ileft,  jleft);
   double s_right = saturation    (iright, jright);
   double c_right = concentration (iright, jright);
   double p_right = pressure      (iright, jright);

   double m_water_left = mobility_water (s_left, c_left);
   double m_oil_left   = mobility_oil (s_left, c_left);
   double m_total_left = m_water_left + m_oil_left;
   double perm_left    = permeability (ileft, jleft);
   double m_perm_left  = m_total_left * perm_left;

   double m_water_right = mobility_water (s_right, c_right);
   double m_oil_right   = mobility_oil (s_right, c_right);
   double m_total_right = m_water_right + m_oil_right;
   double perm_right    = permeability (iright, jright);
   double m_perm_right  = m_total_right * perm_right;

   double m_perm   = harmonic_average (m_perm_left, m_perm_right);

   // we assume dx=dy, so dont divide by length since it is accounted
   // for in flux computation: dont multiply flux by face area
   double dpdn     = (p_right - p_left) / grid.dx;
   double velocity = - m_perm * dpdn;

   // Add gravity effect
   if(g > 0.0)
   {
      double theta_left  = (m_water_left * density_water + 
                            m_oil_left   * density_oil) * gravity * perm_left;
      double theta_right = (m_water_right * density_water + 
                            m_oil_right   * density_oil) * gravity * perm_right;
      double theta = 0.5 * m_perm * ( theta_left/m_perm_left + theta_right/m_perm_right);

      velocity -= theta;
   }

   max_velocity = max ( max_velocity, fabs(velocity) );
   min_velocity = min ( min_velocity, fabs(velocity) );
   
   return velocity;
}

//------------------------------------------------------------------------------
// Find location of minimum for the flux
// This should be called only with g > 0
//------------------------------------------------------------------------------
double argmin_flux(const double& concentration,
                   const double& permeability,
                   const double& velocity)
{
   double r = viscosity_oil / viscosity_water (concentration);
   double z = velocity * viscosity_oil / 
      ((density_water - density_oil) * gravity * permeability);
   double alpha = 27.0 * (-r + r * r - z - 2.0 * r * z - r * r * z);
   double beta  = 2916.0 * r * r * r + alpha * alpha;

   double gamma = pow( alpha + sqrt(beta), 1.0/3.0);
   double fact  = 3.0 * pow(2.0, 1.0/3.0);

   double sstar = (1.0 - fact * r / gamma + gamma / fact) / (1.0 + r);

   double s_min;

   if(sstar <= 1.0 && sstar >= 0.0)
   {
      s_min = sstar;
      ++n_interior_min;
   }
   else
      s_min = (velocity > 0.0) ? 0.0 : 1.0;

   return s_min;
}

//------------------------------------------------------------------------------
// call appropriate numerical flux function
//------------------------------------------------------------------------------
vector<double> ReservoirProblem::num_flux
       (
       const double& velocity,
       const vector<double>& state_left,
       const vector<double>& state_right,
       const double& g
       )
{
   vector<double> flux(2);

   if(flux_type=="dflu")
      flux = dflu_flux (velocity, state_left, state_right, g);

   return flux;
}

//------------------------------------------------------------------------------
// dflu umerical flux function
//------------------------------------------------------------------------------
vector<double> dflu_flux
       (
       const double& velocity,
       const vector<double>& state_left,
       const vector<double>& state_right,
       const double& g
       )
{
   double s_left  = state_left[0];
   double c_left  = state_left[1];
   double s_right = state_right[0];
   double c_right = state_right[1];

   vector<double> flux(2);

   if(g > 0.0) // Gravity is present
   {
      double perm_left = state_left[2];
      double s_min_left = argmin_flux(c_left, perm_left, velocity);

      double perm_right = state_right[2];
      double s_min_right = argmin_flux(c_right, perm_right, velocity);

      s_left  = max( s_left,  s_min_left);
      s_right = min( s_right, s_min_right);

      double m_water_left = mobility_water (s_left, c_left);
      double m_oil_left   = mobility_oil (s_left, c_left);
      double m_total_left = m_water_left + m_oil_left;

      double m_water_right = mobility_water (s_right, c_right);
      double m_oil_right   = mobility_oil (s_right, c_right);
      double m_total_right = m_water_right + m_oil_right;

      double v_left  = velocity 
         - (density_water - density_oil) * gravity * m_oil_left * perm_left;
      double v_right = velocity 
         - (density_water - density_oil) * gravity * m_oil_right * perm_right;
      double f_left  = v_left  * m_water_left / m_total_left;
      double f_right = v_right * m_water_right / m_total_right;
      flux[0]        = max( f_left, f_right );

      if(flux[0] > 0.0)
         flux[1] = c_left  * flux[0];
      else
         flux[1] = c_right * flux[0];
   }
   else // Gravity is not present
   {
      double m_water_left = mobility_water (s_left, c_left);
      double m_oil_left   = mobility_oil (s_left, c_left);
      double m_total_left = m_water_left + m_oil_left;

      double m_water_right = mobility_water (s_right, c_right);
      double m_oil_right   = mobility_oil (s_right, c_right);
      double m_total_right = m_water_right + m_oil_right;

      if (velocity > 0)
      {
         flux[0] = velocity * m_water_left / m_total_left;
         flux[1] = c_left * flux[0];
      }
      else
      {
         flux[0] = velocity * m_water_right / m_total_right;
         flux[1] = c_right * flux[0];
      }
   }

   return flux;
}

//------------------------------------------------------------------------------
// Read some input from file
//------------------------------------------------------------------------------
void ReservoirProblem::read_input ()
{
   cout << "Reading input from file data.in ..." << endl;

   ifstream inp;
   string input;
   inp.open ("data.in");

   inp >> input >> flux_type;      
   assert (input=="flux");
   assert (flux_type=="dflu");

   inp >> input >> order;      
   assert(input == "order");
   assert(order == 1 || order == 2);

   inp >> input >> max_iter;   
   assert(input == "max_iter");
   assert(max_iter > 0);

   inp >> input >> save_freq;   
   assert(input == "save_freq");
   assert(save_freq > 0);

   inp >> input >> cfl;        
   assert(input == "cfl");
   assert(cfl > 0.0 && cfl <= 1.0);

   inp >> input >> cinlet;     
   assert(input == "cinlet");
   assert(cinlet >= 0.0);

   inp >> input >> pinlet;     
   assert(input == "pinlet");

   inp >> input >> poutlet;     
   assert(input == "poutlet");

   inp >> input >> viscosity_oil;     
   assert(input == "mu_oil");
   assert(viscosity_oil > 0.0);

   inp >> input >> density_water;     
   assert(input == "d_water");
   assert(density_water > 0.0);

   inp >> input >> density_oil;     
   assert(input == "d_oil");
   assert(density_oil > 0.0);

   inp >> input >> gravity;
   assert(input == "gravity");
   assert(gravity >= 0.0);

   inp >> input >> grid.xmin >> grid.xmax;
   assert(input == "xrange");

   inp >> input >> grid.ymin >> grid.ymax;
   assert(input == "yrange");

   inp >> grid.nx >> grid.ny;
   inp >> grid.n_boundary;

   // allocate memory for grid
   grid.allocate ();

   for(unsigned int n=0; n<grid.n_boundary; ++n){
      inp >> grid.ibeg[n] >> grid.iend[n]
          >> grid.jbeg[n] >> grid.jend[n]
          >> grid.boundary_condition[n];
      assert (grid.ibeg[n] >= 1 && grid.ibeg[n] <= grid.nx);
      assert (grid.iend[n] >= 1 && grid.iend[n] <= grid.nx);
      assert (grid.jbeg[n] >= 1 && grid.jbeg[n] <= grid.ny);
      assert (grid.jend[n] >= 1 && grid.jend[n] <= grid.ny);
      if(grid.ibeg[n] == grid.iend[n] &&
         grid.jbeg[n] == grid.jend[n])
      {
         cout << "Boundary " << n 
              << " is not a surface !!!" << endl;
         abort ();
      }
   }

   inp.close ();

   cout << "Scheme order           = " << order << endl;
   cout << "Max no. of time steps  = " << max_iter << endl;
   cout << "Solution save freq.    = " << save_freq << endl;
   cout << "CFL number             = " << cfl << endl;
   cout << "Inlet polymer concentr = " << cinlet << endl;
   cout << "Inlet pressure         = " << pinlet << endl;
   cout << "Outlet pressure        = " << poutlet << endl;
   cout << "Viscosity oil          = " << viscosity_oil << endl;
   cout << "Density water          = " << density_water << endl;
   cout << "Density oil            = " << density_oil << endl;
   cout << "Gravity                = " << gravity << endl;
   cout << "nx x ny                = " << grid.nx << " x " << grid.ny << endl;
   cout << "Number of boundaries   = " << grid.n_boundary << endl;
   cout << "Number of cells        = " << grid.n_cells << endl;
   cout << "Number of actual cells = " << (grid.nx-1)*(grid.ny-1) << endl;

}

//------------------------------------------------------------------------------
// Create cartesian grid
//------------------------------------------------------------------------------
void ReservoirProblem::make_grid ()
{
   cout << "Making grid for reservoir problem ..." << endl;

   // set location of boundary
   for(unsigned int n=0; n<grid.n_boundary; ++n)
   {
      if(grid.ibeg[n] == grid.iend[n])
      {
         if(grid.ibeg[n] == 0)
            grid.b_type[n] = imin;
         else
            grid.b_type[n] = imax;
      }

      if(grid.jbeg[n] == grid.jend[n])
      {
         if(grid.jbeg[n] == 0)
            grid.b_type[n] = jmin;
         else
            grid.b_type[n] = jmax;
      }
   }

   grid.dx   = (grid.xmax - grid.xmin)/(grid.nx - 1);
   grid.dy   = (grid.ymax - grid.ymin)/(grid.ny - 1);

   // This is implicitly assumed in the flux computations
   assert (grid.dx == grid.dy);

   // grid vertex coordinates
   for(unsigned int i=0; i<=grid.nx+1; ++i)
      for(unsigned int j=0; j<=grid.ny+1; ++j)
      {
         grid.x (i,j) = grid.xmin + (i-1) * grid.dx;
         grid.y (i,j) = grid.ymin + (j-1) * grid.dy;
      }

   // cell center coordinates
   for(unsigned int i=0; i<=grid.nx; ++i)
      for(unsigned int j=0; j<=grid.ny; ++j)
      {
         grid.xc (i,j) = 0.25 * ( grid.x(i,j)     + grid.x(i+1,j) + 
                                  grid.x(i+1,j+1) + grid.x(i,j+1) );
         grid.yc (i,j) = 0.25 * ( grid.y(i,j)     + grid.y(i+1,j) + 
                                  grid.y(i+1,j+1) + grid.y(i,j+1) );
      }

}

//------------------------------------------------------------------------------
// allocate memory and set initial condition
//------------------------------------------------------------------------------
void ReservoirProblem::initialize ()
{
   saturation.allocate    (grid.nx+1, grid.ny+1);
   concentration.allocate (grid.nx+1, grid.ny+1);
   pressure.allocate      (grid.nx+1, grid.ny+1);
   permeability.allocate  (grid.nx+1, grid.ny+1);

   // rock permeability
   Permeability::N = 50;
   Permeability::xl.resize(Permeability::N);
   Permeability::yl.resize(Permeability::N);
   for(unsigned int i=0; i<Permeability::N; ++i)
   {
      Permeability::xl[i] = (rand() % 1001 + 1) / 1001.0;
      Permeability::yl[i] = (rand() % 1001 + 1) / 1001.0;
   }
   for(unsigned int i=0; i<grid.nx+1; ++i)
      for(unsigned int j=0; j<grid.ny+1; ++j)
         permeability (i,j) = rock_permeability (grid.xc(i,j), grid.yc(i,j));

   // initialize only real cells, not for ghost cells
   for(unsigned int i=1; i<=grid.nx-1; ++i)
      for(unsigned int j=1; j<=grid.ny-1; ++j)
      {
         double dist = grid.xc(i,j) * grid.xc(i,j) + 
                       grid.yc(i,j) * grid.yc(i,j);
         //if(dist <= 0.25 * 0.25) // Water region
         if(grid.yc(i,j) <= -0.25) // Water region
         {
            saturation    (i,j) = 1.0;
            concentration (i,j) = cinlet;
         }
         else // Oil region
         {
            saturation    (i,j) = 0.0;
            concentration (i,j) = 0.0;
         }

         // pressure is everywhere zero to begin with
         pressure (i,j) = 0.0;
      }

   updateGhostCells ();
   output (0);

   final_time = 0.1;

   // RK scheme parameters
   if( order==1 )
      nrk = 1;
   else
      nrk = 3;
   ark[0] = 0.0; ark[1] = 3.0/4.0; ark[2] = 1.0/3.0;
   for (unsigned int i=0; i<nrk; ++i) brk[i] = 1.0 - ark[i];
}

//------------------------------------------------------------------------------
// residual for saturation/concentration equation
//------------------------------------------------------------------------------
void ReservoirProblem::residual (Matrix& s_residual, Matrix& c_residual)
{
   unsigned int i, j;
   double velocity;
   vector<double> state_left(3), state_right(3), flux(2);

   s_residual = 0.0;
   c_residual = 0.0;

   min_velocity = 1.0e20;
   max_velocity = 0.0;

   // interior vertical faces
   for(i=2; i<=grid.nx-1; ++i)
      for(j=1; j<=grid.ny-1; ++j)
      {
         state_left  = reconstruct (i-2, j, i-1, j, i, j);
         state_right = reconstruct (i+1, j, i, j, i-1, j);
         velocity    = darcy_velocity (i-1, j, i, j, 0.0);
         flux        = num_flux (velocity, state_left, state_right, 0.0);

         s_residual (i-1,j) += flux[0] * grid.dy;
         s_residual (i,  j) -= flux[0] * grid.dy;

         c_residual (i-1,j) += flux[1] * grid.dy;
         c_residual (i,  j) -= flux[1] * grid.dy;
      }

   // interior horizontal faces
   for(j=2; j<=grid.ny-1; ++j)
      for(i=1; i<=grid.nx-1; ++i)
      {
         state_left  = reconstruct (i, j-2, i, j-1, i, j);
         state_right = reconstruct (i, j+1, i, j, i, j-1);
         velocity    = darcy_velocity (i, j-1, i, j, gravity);
         flux        = num_flux (velocity, state_left, state_right, gravity);

         s_residual (i,j)   -= flux[0] * grid.dx;
         s_residual (i,j-1) += flux[0] * grid.dx;

         c_residual (i,j)   -= flux[1] * grid.dx;
         c_residual (i,j-1) += flux[1] * grid.dx;
      }

   // inlet/outlet boundaries
   for(unsigned int n=0; n<grid.n_boundary; ++n)
   {
      int bc = grid.boundary_condition[n];

      if (grid.ibeg[n] == grid.iend[n] && bc != SOLID)
      {
         i = grid.ibeg[n];
         for(j=grid.jbeg[n]; j<grid.jend[n]; ++j)
         {

            if (i == 1) // inlet-vertical side
            {
               state_left  = reconstruct (i-1, j, i-1, j, i, j);
               state_right = reconstruct (i+1, j, i, j, i-1, j);
               velocity    = darcy_velocity (i-1, j, i, j, 0.0);
               flux        = num_flux (velocity, state_left, state_right, 0.0);
               s_residual(i,j) -= flux[0] * grid.dy;
               c_residual(i,j) -= flux[1] * grid.dy;
            }
            else // outlet-vertical side
            {
               state_left  = reconstruct (i-2, j, i-1, j, i, j);
               state_right = reconstruct (i, j, i, j, i-1, j);
               velocity    = darcy_velocity (i-1, j, i, j, 0.0);
               flux        = num_flux (velocity, state_left, state_right, 0.0);
               s_residual(i-1,j) += flux[0] * grid.dy;
               c_residual(i-1,j) += flux[1] * grid.dy;
            }
         }
      }

      if (grid.jbeg[n] == grid.jend[n] && bc != SOLID)
      {
         j = grid.jbeg[n];
         for(i=grid.ibeg[n]; i<grid.iend[n]; ++i)
         {

            if(j == 1) // inlet-horizontal side
            {
               state_left  = reconstruct (i, j-1, i, j-1, i, j);
               state_right = reconstruct (i, j+1, i, j, i, j-1);
               velocity    = darcy_velocity (i, j-1, i, j, gravity);
               flux        = num_flux (velocity, state_left, state_right, gravity);
               s_residual(i,j) -= flux[0] * grid.dx;
               c_residual(i,j) -= flux[1] * grid.dx;
            }
            else // outlet-horizontal side
            {
               state_left  = reconstruct (i, j-2, i, j-1, i, j);
               state_right = reconstruct (i, j, i, j, i, j-1);
               velocity    = darcy_velocity (i, j-1, i, j, gravity);
               flux        = num_flux (velocity, state_left, state_right, gravity);
               s_residual(i,j-1) += flux[0] * grid.dx;
               c_residual(i,j-1) += flux[1] * grid.dx;
            }
         }
      }
   }

   dt = cfl * max (grid.dx, grid.dy) / (3.0 * max_velocity);
   double lambda = dt / (grid.dx * grid.dy);
   s_residual *= lambda;
   c_residual *= lambda;
}

//------------------------------------------------------------------------------
// Update polymer concentration
//------------------------------------------------------------------------------
void ReservoirProblem::updateConcentration (Matrix& sc)
{
   for (unsigned int i=1; i<grid.nx; ++i)
      for (unsigned int j=1; j<grid.ny; ++j)
      {
         if (saturation (i,j) > SZERO)
            concentration (i,j) = sc (i,j) / saturation (i,j);
         else
            concentration (i,j) = 0.0;
      }
}

//------------------------------------------------------------------------------
// Update solution in ghost cells
//------------------------------------------------------------------------------
void ReservoirProblem::updateGhostCells ()
{
   unsigned int i, j;

   // top/bottom ghost cells
   for (i=1; i<grid.nx; ++i)
   {
      j = 0;
      saturation    (i,j) = saturation    (i,j+1);
      concentration (i,j) = concentration (i,j+1);
      pressure      (i,j) = pressure      (i,j+1);

      j = grid.ny;
      saturation    (i,j) = saturation    (i,j-1);
      concentration (i,j) = concentration (i,j-1);
      pressure      (i,j) = pressure      (i,j-1);
   }

   // left/right ghost cells
   for (j=1; j<grid.ny; ++j)
   {
      i = 0;
      saturation    (i,j) = saturation    (i+1,j);
      concentration (i,j) = concentration (i+1,j);
      pressure      (i,j) = pressure      (i+1,j);

      i = grid.nx;
      saturation    (i,j) = saturation    (i-1,j);
      concentration (i,j) = concentration (i-1,j);
      pressure      (i,j) = pressure      (i-1,j);
   }

   // inlet/outlet boundaries
   for(unsigned int n=0; n<grid.n_boundary; ++n)
   {
      if (grid.ibeg[n] == grid.iend[n])
      {
         i = grid.ibeg[n];
         for(j=grid.jbeg[n]; j<grid.jend[n]; ++j)
         {

            if (grid.ibeg[n] == 1) // inlet-vertical side
            {
               saturation    (i-1,j) = 1.0;
               concentration (i-1,j) = cinlet;
               pressure      (i-1,j) = pinlet;
            }
            else // outlet-vertical side
            {
               //saturation    (i,j) = 0.0;
               //concentration (i,j) = 0.0;
               pressure      (i,j) = poutlet;
            }
         }
      }

      if (grid.jbeg[n] == grid.jend[n])
      {
         j = grid.jbeg[n];
         for(i=grid.ibeg[n]; i<grid.iend[n]; ++i)
         {

            if(grid.jbeg[n] == 1) // inlet-horizontal side
            {
               saturation    (i,j-1) = 1.0;
               concentration (i,j-1) = cinlet;
               pressure      (i,j-1) = pinlet;
            }
            else // outlet-horizontal side
            {
               //saturation    (i,j) = 1.0;
               //concentration (i,j) = 0.0;
               pressure      (i,j) = poutlet;
            }
         }
      }
   }

}

//------------------------------------------------------------------------------
// Find min and max values of solution
//------------------------------------------------------------------------------
void ReservoirProblem::findMinMax () const
{
   double s_min = 1.0e20;
   double s_max =-1.0e20;
   double c_min = 1.0e20;
   double c_max =-1.0e20;
   double p_min = 1.0e20;
   double p_max =-1.0e20;

   for (unsigned int i=1; i<grid.nx; ++i)
      for (unsigned int j=1; j<grid.ny; ++j)
      {
         s_min = min (s_min, saturation(i,j));
         s_max = max (s_max, saturation(i,j));
         c_min = min (c_min, concentration(i,j));
         c_max = max (c_max, concentration(i,j));
         p_min = min (p_min, pressure(i,j));
         p_max = max (p_max, pressure(i,j));
      }

   cout << "Saturation    = " << s_min << " " << s_max << endl;
   cout << "Concentration = " << c_min << " " << c_max << endl;
   cout << "Pressure      = " << p_min << " " << p_max << endl;
   cout << "Min velocity  = " << min_velocity << endl;
   cout << "Max velocity  = " << max_velocity << endl;
   cout << "dt            = " << dt << endl;
}

//------------------------------------------------------------------------------
// perform time stepping
//------------------------------------------------------------------------------
void ReservoirProblem::solve ()
{
   unsigned int iter = 0;
   double time = 0.0;
   PressureProblem pressure_problem (&grid);
   Matrix s_residual (grid.nx+1, grid.ny+1);
   Matrix c_residual (grid.nx+1, grid.ny+1);
   Matrix sc         (grid.nx+1, grid.ny+1);
   Matrix s_old      (grid.nx+1, grid.ny+1);
   Matrix sc_old     (grid.nx+1, grid.ny+1);

   while (iter < max_iter)
   {
      s_old = saturation;
      sc_old= saturation * concentration;

      // solve for pressure
      pressure_problem.run (saturation, concentration, 
                            permeability, pressure);

      // Runge-Kutta stages
      for (unsigned int irk=0; irk<nrk; ++irk)
      {
         n_interior_min = 0; // reset counter

         // compute residual
         residual (s_residual, c_residual);
      
         // new value of s*c
         sc = sc_old * ark[irk] + 
              (saturation * concentration - c_residual) * brk[irk];

         // update saturation
         saturation  = s_old * ark[irk] + (saturation - s_residual) * brk[irk];

         // update concentration
         updateConcentration (sc);

         // update solution in ghost cells
         updateGhostCells ();
      }

      // find solution range: to check for stability
      findMinMax ();
      
      time += dt;
      ++iter;

      // save solution to file
      if (iter % save_freq == 0 || iter == max_iter)
         output (iter);

      cout << "Time= " << time << " iter= " << iter << endl;
      cout << "No. of interior min flux = " << n_interior_min << endl;
      cout << endl;
   }
}

//------------------------------------------------------------------------------
// save solution to file
// only interior cells are written, ghost cells are not written
//------------------------------------------------------------------------------
void ReservoirProblem::output (const unsigned int iter) const
{

   unsigned int i, j;
   ofstream vtk;
   ostringstream filename;
   filename << "solution-" << iter << ".vtk";

   vtk.open (filename.str().c_str());

   vtk << "# vtk DataFile Version 2.0" << endl;
   vtk << "Oil Reservoir Problem: iter = " << iter << endl;
   vtk << "ASCII" << endl;
   vtk << "DATASET STRUCTURED_GRID" << endl;
   vtk << "DIMENSIONS " << grid.nx << " " << grid.ny << " 1" << endl;

   // write coordinates
   vtk << "POINTS " << grid.nx * grid.ny << " float" << endl;
   for(j=1; j<=grid.ny; ++j)
      for(i=1; i<=grid.nx; ++i)
      {
         vtk << grid.x (i,j) << "  "
             << grid.y (i,j) << "  "
             << 0.0 << endl;
      }

   vtk << "CELL_DATA " << (grid.nx-1)*(grid.ny-1) << endl;

   // write saturation
   vtk << "SCALARS saturation float 1" << endl;
   vtk << "LOOKUP_TABLE default" << endl;
   for(j=1; j<grid.ny; ++j)
      for(i=1; i<grid.nx; ++i)
         vtk << fixed << saturation (i,j) << endl;

   // write concentration
   vtk << "SCALARS concentration float 1" << endl;
   vtk << "LOOKUP_TABLE default" << endl;
   for(j=1; j<grid.ny; ++j)
      for(i=1; i<grid.nx; ++i)
         vtk << fixed << concentration (i,j) << endl;

   // write pressure
   vtk << "SCALARS pressure float 1" << endl;
   vtk << "LOOKUP_TABLE default" << endl;
   for(j=1; j<grid.ny; ++j)
      for(i=1; i<grid.nx; ++i)
         vtk << fixed << pressure (i,j) << endl;

   // write permeability
   if (iter==0)
   {
      vtk << "SCALARS permeability float 1" << endl;
      vtk << "LOOKUP_TABLE default" << endl;
      for(j=1; j<grid.ny; ++j)
         for(i=1; i<grid.nx; ++i)
            vtk << fixed << permeability (i,j) << endl;
   }

   vtk.close ();

}

//------------------------------------------------------------------------------
// solve the whole problem
//------------------------------------------------------------------------------
void ReservoirProblem::run ()
{
   read_input ();
   make_grid ();
   initialize ();
   solve ();
}
