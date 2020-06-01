#include <iostream>
#include <valarray>
#include "matrix.h"
#include "grid.h"
#include "pressure.h"
#include "material.h"

using namespace std;

//------------------------------------------------------------------------------
// constructor given grid
//------------------------------------------------------------------------------
PressureProblem::PressureProblem (Grid* grid_in)
{
   // set pointer grid to grid_in
   grid = grid_in;
}

//------------------------------------------------------------------------------
// compute rhs (b) of pressure equation A*p=b
//------------------------------------------------------------------------------
Matrix PressureProblem::compute_rhs (const Matrix& saturation,
                                     const Matrix& concentration,
                                     const Matrix& permeability,
                                     const Matrix& pressure)
{
   unsigned int i, j;
   double mobility_water_left, mobility_oil_left;
   double mobility_water_right, mobility_oil_right;
   double mobility_left, mobility_right;
   double perm_left, perm_right;
   double theta_left, theta_right, theta;
   double m_perm_left, m_perm_right, m_perm;
   double flux;
   Matrix result(grid->nx+1, grid->ny+1);

   result = 0.0;

   // interior horizontal faces
   // contribution from gravity term
   for(j=2; j<=grid->ny-1; ++j)
      for(i=1; i<=grid->nx-1; ++i)
      {
         mobility_water_left = mobility_water (saturation(i,j-1), concentration(i,j-1));
         mobility_oil_left = mobility_oil (saturation(i,j-1), concentration(i,j-1));
         mobility_left = mobility_water_left + mobility_oil_left;
         perm_left = permeability (i,j-1);
         m_perm_left = mobility_left * perm_left;
         theta_left = (mobility_water_left * density_water +
                       mobility_oil_left   * density_oil) * gravity * perm_left;

         mobility_water_right = mobility_water (saturation(i,j), concentration(i,j));
         mobility_oil_right = mobility_oil (saturation(i,j), concentration(i,j));
         mobility_right = mobility_water_right + mobility_oil_right;
         perm_right = permeability (i,j);
         m_perm_right = mobility_right * perm_right;
         theta_right = (mobility_water_right * density_water +
                        mobility_oil_right   * density_oil) * gravity * perm_right;

         m_perm = harmonic_average (m_perm_left, m_perm_right);

         theta  = 0.5 * m_perm * ( theta_left/m_perm_left + theta_right/m_perm_right );
         flux   = theta * grid->dx;

         result(i,j)   -= flux;
         result(i,j-1) += flux;
      }

   // inlet/outlet boundaries
   for(unsigned int n=0; n<grid->n_boundary; ++n)
   {
      int bc = grid->boundary_condition[n];

      if (grid->ibeg[n] == grid->iend[n]) // Vertical boundary faces
      {
         i = grid->ibeg[n];
         for(j=grid->jbeg[n]; j<grid->jend[n]; ++j)
         {
            mobility_left = mobility_total (saturation(i-1,j), concentration(i-1,j));
            perm_left = permeability (i-1,j);
            mobility_right = mobility_total (saturation(i,j), concentration(i,j));
            perm_right = permeability (i,j);
            m_perm = harmonic_average (mobility_left  * perm_left, 
                                       mobility_right * perm_right);

            if (i == 1 && bc == INLET) // inlet-vertical side
            {
               // dpdn = (pressure(i,j) - pinlet)/(dx)
               flux         = m_perm * (-pinlet)/(grid->dx) * grid->dy;
               result(i,j) -= flux;
            }
            else if (i == grid->nx && bc == OUTLET) // outlet-vertical side
            {
               // dpdn = (poutlet - pressure(i-1,j))/(dx)
               flux           = m_perm * (poutlet)/(grid->dx) * grid->dy;
               result(i-1,j) += flux;
            }
            else
            {
               cout << "pressure rhs: Boundary index is wrong !!!\n";
               abort ();
            }
         }
      }

      if (grid->jbeg[n] == grid->jend[n]) // Horizontal boundary faces
      {
         j = grid->jbeg[n];
         for(i=grid->ibeg[n]; i<grid->iend[n]; ++i)
         {
            mobility_water_left = mobility_water (saturation(i,j-1), concentration(i,j-1));
            mobility_oil_left = mobility_oil (saturation(i,j-1), concentration(i,j-1));
            mobility_left = mobility_water_left + mobility_oil_left;
            perm_left = permeability (i,j-1);
            m_perm_left = mobility_left * perm_left;
            theta_left = (mobility_water_left * density_water +
                          mobility_oil_left   * density_oil) * gravity * perm_left;

            mobility_water_right = mobility_water (saturation(i,j), concentration(i,j));
            mobility_oil_right = mobility_oil (saturation(i,j), concentration(i,j));
            mobility_right = mobility_water_right + mobility_oil_right;
            perm_right = permeability (i,j);
            m_perm_right = mobility_right * perm_right;
            theta_right = (mobility_water_right * density_water +
                           mobility_oil_right   * density_oil) * gravity * perm_right;

            m_perm = harmonic_average (m_perm_left, m_perm_right);

            theta  = 0.5 * m_perm * ( theta_left/m_perm_left + theta_right/m_perm_right );

            if(j == 1 && bc == INLET) // inlet-horizontal side
            {
               // dpdn = (pressure(i,j) - pinlet)/(dy)
               flux         = m_perm * (-pinlet)/(grid->dy) * grid->dx
                            + theta * grid->dx;
               result(i,j) -= flux;
            }
            else if(j == grid->ny && bc == OUTLET) // outlet-horizontal side
            {
               // dpdn = (poutlet - pressure(i,j-1))/(dy)
               flux           = m_perm * (poutlet)/(grid->dy) * grid->dx
                              + theta * grid->dx;
               result(i,j-1) += flux;
            }
            else
            {
               cout << "pressure rhs: Boundary index is wrong !!!\n";
               abort ();
            }
         }
      }
   }

   return result;

}

//------------------------------------------------------------------------------
// compute matrix vector product A*p in pressure equation A*p = b
//------------------------------------------------------------------------------
Matrix PressureProblem::A_times_pressure (const Matrix& saturation,
                                          const Matrix& concentration,
                                          const Matrix& permeability,
                                          const Matrix& pressure)
{
   unsigned int i, j;
   double mobility_left, mobility_right;
   double perm_left, perm_right;
   double m_perm, dpdn, flux;
   Matrix result(grid->nx+1, grid->ny+1);

   result = 0.0;

   // interior vertical faces
   for(i=2; i<=grid->nx-1; ++i)
      for(j=1; j<=grid->ny-1; ++j)
      {
         mobility_left = mobility_total (saturation(i-1,j), concentration(i-1,j));
         perm_left = permeability (i-1,j);
         mobility_right = mobility_total (saturation(i,j), concentration(i,j));
         perm_right = permeability (i,j);
         m_perm = harmonic_average (mobility_left  * perm_left, 
                                    mobility_right * perm_right);

         dpdn = (pressure(i,j) - pressure(i-1,j))/grid->dx;

         flux           = - m_perm * dpdn * grid->dy;
         result(i-1,j) += flux;
         result(i,j)   -= flux;
      }

   // interior horizontal faces
   for(j=2; j<=grid->ny-1; ++j)
      for(i=1; i<=grid->nx-1; ++i)
      {
         mobility_left = mobility_total (saturation(i,j-1), concentration(i,j-1));
         perm_left = permeability (i,j-1);
         mobility_right = mobility_total (saturation(i,j), concentration(i,j));
         perm_right = permeability (i,j);
         m_perm = harmonic_average (mobility_left  * perm_left, 
                                    mobility_right * perm_right);

         dpdn = (pressure(i,j) - pressure(i,j-1))/grid->dy;

         flux           = - m_perm * dpdn * grid->dx;
         result(i,j)   -= flux;
         result(i,j-1) += flux;
      }

   // inlet/outlet boundaries: vertical
   for(unsigned int n=0; n<grid->n_boundary; ++n)
   {
      if (grid->ibeg[n] == grid->iend[n])
      {
         i = grid->ibeg[n];
         for(j=grid->jbeg[n]; j<grid->jend[n]; ++j)
         {
            mobility_left = mobility_total (saturation(i-1,j), concentration(i-1,j));
            perm_left = permeability (i-1,j);
            mobility_right = mobility_total (saturation(i,j), concentration(i,j));
            perm_right = permeability (i,j);
            m_perm = harmonic_average (mobility_left  * perm_left, 
                                       mobility_right * perm_right);

            if (grid->ibeg[n] == 1) // inlet-vertical side
            {
               // dpdn = (pressure(i,j) - pinlet)/(dx)
               flux         = - m_perm * pressure(i,j)/(grid->dx) * grid->dy;
               result(i,j) -= flux;
            }
            else // outlet-vertical side
            {
               // dpdn = (poutlet - pressure(i-1,j))/(dx)
               flux           = - m_perm * (-pressure(i-1,j))/(grid->dx) * grid->dy;
               result(i-1,j) += flux;
            }
         }
      }

      // inlet/outlet boundaries: horizontal
      if (grid->jbeg[n] == grid->jend[n])
      {
         j = grid->jbeg[n];
         for(i=grid->ibeg[n]; i<grid->iend[n]; ++i)
         {
            mobility_left = mobility_total (saturation(i,j-1), concentration(i,j-1));
            perm_left = permeability (i,j-1);
            mobility_right = mobility_total (saturation(i,j), concentration(i,j));
            perm_right = permeability (i,j);
            m_perm = harmonic_average (mobility_left  * perm_left, 
                                       mobility_right * perm_right);

            if(grid->jbeg[n] == 1) // inlet-horizontal side
            {
               // dpdn = (pressure(i,j) - pinlet)/(dy)
               flux         = - m_perm * (pressure(i,j))/(grid->dy) * grid->dx;
               result(i,j) -= flux;
            }
            else // outlet-horizontal side
            {
               // dpdn = (poutlet - pressure(i,j-1))/(dy)
               flux           = - m_perm * (-pressure(i,j-1))/(grid->dy) * grid->dx;
               result(i,j-1) += flux;
            }
         }
      }
   }

   return result;

}

//------------------------------------------------------------------------------
// Compute residual for pressure problem, r = b - A*p
//------------------------------------------------------------------------------
Matrix PressureProblem::residual (const Matrix& saturation, 
                                  const Matrix& concentration,
                                  const Matrix& permeability,
                                  const Matrix& pressure)
{
   Matrix r(grid->nx+1, grid->ny+1);

   r = compute_rhs (saturation, concentration, permeability, pressure)
     - A_times_pressure(saturation, concentration, permeability, pressure);

   return r;
}

//------------------------------------------------------------------------------
// Solve pressure equation by CG method
//------------------------------------------------------------------------------
void PressureProblem::run (const Matrix& saturation, 
                           const Matrix& concentration,
                           const Matrix& permeability,
                                 Matrix& pressure)
{
   const unsigned int max_iter = 5000;
   const double tolerance = 1.0e-6;
   unsigned int iter = 0;
   double beta, omega;
   vector<double> r2(max_iter);
   Matrix d (grid->nx + 1, grid->ny + 1);
   Matrix r (grid->nx + 1, grid->ny + 1);
   Matrix v (grid->nx + 1, grid->ny + 1);

   // initial residual
   r = residual (saturation, concentration, permeability, pressure);

   // initial direction
   d = r;

   r2[0] = r.dot(r);

   // CG iterations
   while ( sqrt(r2[iter]/r2[0]) > tolerance && iter < max_iter )
   {
      if (iter >= 1) // update descent direction
      {              // d = r + beta * d
         beta = r2[iter] / r2[iter-1];
         d   *= beta;
         d   += r;
      }

      v = A_times_pressure (saturation, concentration, permeability, d);
      omega = r2[iter] / d.dot(v);

      // update pressure: p = p + omega * d
      pressure += d * omega;

      // update residual: r = r - omega * v
      v *= omega;
      r -= v;

      ++iter;

      r2[iter] = r.dot(r);
   }

   cout << "PressureProblem: iter= " << iter 
        << " residue= " << sqrt(r2[iter]/r2[0]) << endl;

   if (sqrt(r2[iter]) > tolerance && iter==max_iter)
   {
      cout << "PressureProblem did not converge !!!" << endl;
      abort ();
   }
}
