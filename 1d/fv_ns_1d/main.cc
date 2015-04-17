#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#define SIGN(a) (((a)<0) ? -1:1)

const double Pr    = 2.0/3.0;
const double arks[] = {0.0, 3.0/4.0, 1.0/3.0};
const double brks[] = {1.0, 1.0/4.0, 2.0/3.0};

// These values are set below based on test case type
double GAMMA;
double gas_const;
double mu_ref, T_ref;

using namespace std;

// viscosity coefficient as functon of temperature
double viscosity (const double T)
{
   double mu = mu_ref * pow( T / T_ref, 0.8);
   return mu;
}

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
vector<double> muscl (const vector<double>& ul,
                      const vector<double>& uc,
                      const vector<double>& ur)
{
   unsigned int n = ul.size();
   vector<double> result (n);
   double dul, duc, dur;
   const double beta = 2.0;
   
   for(unsigned int i=0; i<n; ++i)
   {
      dul = uc[i] - ul[i];
      dur = ur[i] - uc[i];
      duc = (ur[i] - ul[i])/2.0;
      //result[i] = uc[i] + 0.5 * vanleer (dul, dur);
      result[i] = uc[i] + 0.5 * minmod (beta*dul, duc, beta*dur);
   }
   
   return result;
}

//------------------------------------------------------------------------------
// Main class of the problem
//------------------------------------------------------------------------------
class FVProblem
{
   public:
      FVProblem ();
      void run ();
   
   private:
      void make_grid_and_dofs ();
      void initial_condition ();
      void compute_dt ();
      void con_to_prim ();
      void compute_face_derivatives ();
      void reconstruct (const unsigned int face,
                        vector<double>& left,
                        vector<double>& right) const;
      void kfvs_split_flux (const vector<double>& prim,
                            const double& tau,
                            const double& q,
                            const int sign,
                            vector<double>& flux) const;
      void split_U (const vector<double>& prim,
                    const int sign,
                    vector<double>& U) const;
      void rkfvs (const vector<double>&,
                  const vector<double>&,
                        vector<double>&) const;
      void num_flux (const vector<double>&,
                     const vector<double>&,
                     const double&         tau_left,
                     const double&         tau_right,
                     const double&         q_left,
                     const double&         q_right,
                     vector<double>&) const;
      void compute_residual ();
      void update_solution (const unsigned int rk);
      void output ();
      
      double d_left, u_left, p_left;
      double d_right, u_right, p_right;
      vector<double> prim_left;
      vector<double> prim_right;
      unsigned int n_var;
      unsigned int n_cell;
      unsigned int n_face;
      double xmin;
      double xmax;
      double xmid;
      double dx;
      double dt;
      double cfl;
      double final_time;
      vector<double> xc;
      vector<double> xf;
   
      vector< vector<double> > primitive;
      vector< vector<double> > residual;
      vector< vector<double> > conserved;
      vector< vector<double> > conserved_old;
      vector<double> tau_f, q_f;

};

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
FVProblem::FVProblem ()
{   
   n_var  = 3;

   int test_case = 1;

   if(test_case == 1)
   {
      // Sod shock tube case
      GAMMA = 1.4;
      gas_const = 1.0;
      final_time = 0.2;
      n_cell = 100;
      cfl    = 0.8;

      xmin    = 0.0;
      xmax    = 1.0;
      xmid   = 0.5 * (xmin + xmax);

      d_left  = 1.0;
      d_right = 0.125;
   
      u_left  = 0.0;
      u_right = 0.0;
   
      p_left  = 1.0;
      p_right = 0.1;

      mu_ref= 0.0;
      T_ref = p_left / d_left / gas_const;
   }
   else if(test_case == 2)
   {
      // shock structure case
      GAMMA = 5.0/3.0;
      gas_const = 0.5;
      final_time = 200.0;
      n_cell = 100;
      cfl    = 0.1;

      xmin   = -0.25;
      xmax   =  0.00;
      xmid   = 0.5 * (xmin + xmax);
      double mach_left = 1.5;
      double M2 = pow(mach_left, 2);

      d_left = 1.0;
      u_left = 1.0;
      p_left = 1.0/GAMMA/M2;

      d_right= (GAMMA+1.0)*M2/(2.0+(GAMMA-1.0)*M2)*d_left;
      u_right= ((GAMMA-1.0)/(GAMMA+1.0)+2.0/(GAMMA+1.0)/M2)*u_left;
      p_right= (2.0*GAMMA/(GAMMA+1.0)*M2-(GAMMA-1)/(GAMMA+1.0))*p_left;

      mu_ref= 0.0005;
      T_ref = p_left / d_left / gas_const;
   }

   prim_left.resize (n_var);
   prim_right.resize (n_var);

   prim_left[0] = d_left;
   prim_left[1] = u_left;
   prim_left[2] = p_left;
   
   prim_right[0] = d_right;
   prim_right[1] = u_right;
   prim_right[2] = p_right;
   
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
   
   primitive.resize (n_cell, vector<double>(n_var));
   residual.resize (n_cell, vector<double>(n_var));
   conserved.resize (n_cell, vector<double>(n_var));
   conserved_old.resize (n_cell, vector<double>(n_var));
   
   tau_f.resize (n_face);
   q_f.resize   (n_face);
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
      {
         conserved[i][0] = d_left;
         conserved[i][1] = d_left * u_left;
         conserved[i][2] = p_left/(GAMMA-1.0) + 0.5 * d_left * pow(u_left,2);
      }
      else
      {
         conserved[i][0] = d_right;
         conserved[i][1] = d_right * u_right;
         conserved[i][2] = p_right/(GAMMA-1.0) + 0.5 * d_right * pow(u_right,2);
      }
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
      double speed = fabs(primitive[i][1]) + 
                     sqrt(GAMMA * primitive[i][2] / primitive[i][0]);
      dt = min (dt, dx/speed);
   }

   dt = min (dt, dx * dx / mu_ref);
   dt *= cfl;
}

//------------------------------------------------------------------------------
// Convert conserved to primitive
//------------------------------------------------------------------------------
void FVProblem::con_to_prim ()
{
   for(unsigned int i=0; i<n_cell; ++i)
   {
      primitive[i][0] = conserved[i][0];
      primitive[i][1] = conserved[i][1]/conserved[i][0];
      primitive[i][2] = (GAMMA-1.0) * (conserved[i][2] - 
                           0.5 * pow(conserved[i][1], 2.0) / conserved[i][0]);
   }
}

//------------------------------------------------------------------------------
// Compute derivatives at faces
//------------------------------------------------------------------------------
void FVProblem::compute_face_derivatives ()
{
   tau_f[0] = 0.0;
   q_f[0]   = 0.0;
   
   for(unsigned int i=1; i<n_face-1; ++i)
   {
      double T_left  = primitive[i-1][2] / (primitive[i-1][0] * gas_const);
      double T_right = primitive[i][2] / (primitive[i][0] * gas_const);
      double T       = 0.5 * (T_left + T_right);
      double mu      = viscosity (T);

      tau_f[i] = (4.0 * mu / 3.0 ) * 
                 (primitive[i][1] - primitive[i-1][1]) / dx;

      double k = mu * GAMMA * gas_const / (GAMMA-1.0) / Pr;

      q_f[i]   = -k * (T_right - T_left) / dx;
   }
   
   tau_f[n_face-1] = 0.0;
   q_f[n_face-1]   = 0.0;
}

//------------------------------------------------------------------------------
// Reconstruct left/right state at a face
//------------------------------------------------------------------------------
void FVProblem::reconstruct (const unsigned int face,
                             vector<double>& left,
                             vector<double>& right) const
{
   if(face==1)
   {
      left = primitive[0];
      right = muscl (primitive[2], primitive[1], primitive[0]);
   }
   else if(face==n_face-2)
   {
      left  = muscl (primitive[n_cell-3], primitive[n_cell-2], 
                     primitive[n_cell-1]);
      right = primitive[n_cell-1];
   }
   else
   {
      left  = muscl (primitive[face-2], primitive[face-1], primitive[face]);
      right = muscl (primitive[face+1], primitive[face], primitive[face-1]);
   }
}

//------------------------------------------------------------------------------
// KFVS split fluxes: sign=+1 give positive flux and
// sign=-1 gives negative flux
//------------------------------------------------------------------------------
void FVProblem::kfvs_split_flux (const vector<double>& prim,
                                 const double& tau,
                                 const double& q,
                                 const int sign,
                                 vector<double>& flux) const
{
   double beta, s, A, B, E, fact;

   beta = 0.5 * prim[0] / prim[2];
   s    = prim[1] * sqrt(beta);
   A    = 0.5 * (1.0 + sign * erf(s));
   B    = sign * 0.5 * exp(-s * s) / sqrt(beta * M_PI);
   E    = prim[2]/(GAMMA-1.0) + 0.5 * prim[0] * pow(prim[1], 2);
   fact = prim[1] * A + B;
   
   // inviscid flux
   flux[0] = prim[0] * fact;
   flux[1] = (prim[2] + prim[0] * pow(prim[1], 2)) * A + 
             prim[0] * prim[1] * B;
   flux[2] = prim[1] * (E + prim[2]) * A +
             (E + 0.5 * prim[2]) * B;

   double g1 = (GAMMA - 1.0) / GAMMA;
   double g2 = (3.0 * GAMMA - 1.0) / (4.0 * (GAMMA - 1.0) );
   double g3 = (3.0 - GAMMA) / (2.0 * GAMMA);

   // viscous flux
   flux[0] += prim[0] * B / prim[2] * ( -0.5 * tau -
                  g1 * prim[1] * beta * q );
   flux[1] += -tau * A + g1 * prim[0] * B * q / prim[2];
   flux[2] += (-prim[1] * tau + q) * A +
              B * (-g2 * tau - g3 * prim[1] * beta * q);
}

//------------------------------------------------------------------------------
// KFVS split conserved variables: 
//------------------------------------------------------------------------------
void FVProblem::split_U (const vector<double>& prim,
                         const int sign,
                         vector<double>& U) const
{
   double beta, s, A, B, E, fact;

   beta = 0.5 * prim[0] / prim[2];
   s    = prim[1] * sqrt(beta);
   A    = 0.5 * (1.0 + sign * erf(s));
   B    = sign * 0.5 * exp(-s * s) / sqrt(beta * M_PI);
   E    = prim[2]/(GAMMA-1.0) + 0.5 * prim[0] * pow(prim[1], 2);
   fact = prim[1] * A + B;
   
   // inviscid flux
   U[0] = prim[0] * A;
   U[1] = prim[0] * fact;
   U[2] = E * A + 0.5 * prim[0] * prim[1] * B;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void FVProblem::rkfvs(const vector<double>& left,
                      const vector<double>& right,
                            vector<double>& flux) const
{
   vector<double> Up(n_var);
   vector<double> Um(n_var);
   vector<double> U(n_var);

   split_U(left,  +1,  Up);
   split_U(right, -1, Um);

   for(unsigned int i=0; i<n_var; ++i)
      U[i] = Up[i] + Um[i];

   double rho = U[0];
   double u   = U[1] / U[0];
   double p   = (GAMMA-1.0) * (U[2] - 0.5 * rho * u * u);

   flux[0] = rho * u;
   flux[1] = p + rho * u * u;
   flux[2] = (U[2] + p) * u;
}

//------------------------------------------------------------------------------
// Numerical flux function
//------------------------------------------------------------------------------
void FVProblem::num_flux(const vector<double>& left,
                         const vector<double>& right,
                         const double&         tau_left,
                         const double&         tau_right,
                         const double&         q_left,
                         const double&         q_right,
                         vector<double>&       flux) const
{
   //rkfvs(left, right, flux); return;

   vector<double> flux_pos(n_var);
   vector<double> flux_neg(n_var);
   
   kfvs_split_flux (left,  tau_left,  q_left,  +1, flux_pos);
   kfvs_split_flux (right, tau_right, q_right, -1, flux_neg);
   
   for(unsigned int i=0; i<n_var; ++i)
      flux[i] = flux_pos[i] + flux_neg[i];
   
  
}

//------------------------------------------------------------------------------
// Compute finite volume residual
//------------------------------------------------------------------------------
void FVProblem::compute_residual ()
{
   for(unsigned int i=0; i<n_cell; ++i)
      for(unsigned int j=0; j<n_var; ++j)
         residual[i][j] = 0.0;
   
   vector<double> flux (n_var);
   vector<double> left (n_var);
   vector<double> right(n_var);
   
   // Flux through left boundary
   num_flux (prim_left, primitive[0], 
             tau_f[0], tau_f[0], 
             q_f[0], q_f[0], flux);
   for(unsigned int j=0; j<n_var; ++j)
      residual[0][j] -= flux[j];

   // Flux through interior faces
   for(unsigned int i=1; i<n_face-1; ++i)
   {
      reconstruct (i, left, right);
      num_flux (left, right, tau_f[i], tau_f[i], q_f[i], q_f[i], flux);
      for(unsigned int j=0; j<n_var; ++j)
      {
         residual[i-1][j] += flux[j];
         residual[i][j]   -= flux[j];
      }
   }

   // Flux through right boundary
   num_flux (primitive[n_cell-1], prim_right, 
             tau_f[n_face-1], tau_f[n_face-1], 
             q_f[n_face-1], q_f[n_face-1], flux);
   for(unsigned int j=0; j<n_var; ++j)
      residual[n_cell-1][j] += flux[j];
}

//------------------------------------------------------------------------------
// Compute finite volume residual
//------------------------------------------------------------------------------
void FVProblem::update_solution (const unsigned int rk)
{
   for(unsigned int i=0; i<n_cell; ++i)
      for(unsigned int j=0; j<n_var; ++j)
         conserved[i][j] = arks[rk] * conserved_old[i][j] +
            brks[rk] * (conserved[i][j] - (dt/dx) * residual[i][j]);
   
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
         << primitive[i][0] << " " 
         << primitive[i][1] << " "
         << primitive[i][2] << endl;
   fo.close ();
}

//------------------------------------------------------------------------------
// Start the computations
//------------------------------------------------------------------------------
void FVProblem::run ()
{
   make_grid_and_dofs ();
   initial_condition ();
   con_to_prim ();

   double time = 0.0;
   unsigned int iter = 0;
   while (time < final_time)
   {
      conserved_old = conserved;
      compute_dt ();
      if(time+dt > final_time) dt = final_time - time;
      for(unsigned int rk=0; rk<3; ++rk)
      {
         compute_face_derivatives ();
         compute_residual ();
         update_solution (rk);
         con_to_prim ();
      }
      time += dt;
      ++iter;
      if(iter % 1000 == 0) 
      cout << "Iter = " << iter << " Time = " << time << endl;
      if(iter % 1000 == 0) output ();
   }
   
   con_to_prim ();
   output ();
}

//------------------------------------------------------------------------------
int main ()
{
   FVProblem fv_problem;
   fv_problem.run ();
   
   return 0;
}
