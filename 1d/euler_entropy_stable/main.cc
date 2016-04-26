#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

enum ReconstructionScheme { FIRST, MINMOD, VANLEER, WENO};

enum FluxScheme {KEPSSENT, ROE, ROEFIXED, RUSANOV};

enum TimeIntegrationScheme{RK1, SSPRK3, JAMESON_RK4};

#define SIGN(a) (((a)<0) ? -1:1)
#define Cp  (GAMMA * gas_const / (GAMMA - 1.0))

// Coefficients for RK scheme
const double arks[] = {0.0, 3.0/4.0, 1.0/3.0};
const double brks[] = {1.0, 1.0/4.0, 2.0/3.0};
const double jameson_rks[] = {1.0/4.0,1.0/3.0,1.0/2.0, 1.0};

// These values are set below based on test case type
double GAMMA;
double gas_const;
double K2, K4;
double beta_upwind,alpha,beta_upwind_cen,alpha_cen; // factor for increasing wave speed
double nrk; // no. of rk stages: 1 or 3 only

unsigned int counter = 0; // For solution storage

using namespace std;

void con2prim (vector<double>& con, vector<double>& prim);



//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
double logavg(double a, double b)
{
   double xi = b/a;
   double f = (xi - 1.0) / (xi + 1.0);
   double u = f * f;
   
   double F;
   if (u < 1.0e-6)
   {
      double u2 = u * u;
      double u3 = u2 * u;
      F = 1.0 + u/3.0 + u2/5.0 + u3/7.0;
   }
   else
      F = log(xi)/2.0/f;
   
   return 0.5*(a+b)/F;
}


//------------------------------------------------------------------------------
// Compute temperature given primitive variables
//------------------------------------------------------------------------------
double temperature(const vector<double>& prim)
{
   return prim[2] / (gas_const * prim[0]);
}

//------------------------------------------------------------------------------
// Compute temperature given primitive variables
//------------------------------------------------------------------------------
double enthalpy(const vector<double>& prim)
{
   return GAMMA * prim[2] / (prim[0] * (GAMMA-1.0)) + 
          0.5 * pow(prim[1], 2);
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
// Minmod for entropy variable (treated like an indicator function)
//------------------------------------------------------------------------------
double minmod2 (const double& a,
               const double& b)
{
   double result;
   
   if (a*b < 0.0)
      result = 0.0;
   else if (a < 0.0 && b == 0.0)
      result = 0.0;
   else if (a >= 0.0 && b == 0.0)
      result = 0.0;   
   else if (a/b >= 1)
      result = 1.0;   
   else        
      result = a/b;
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
// Weno reconstruction
//------------------------------------------------------------------------------
double weno5(double um2, double um1, double u0, double up1, double up2)
{
   double eps = 1.0e-6;
   double gamma1=1.0/10.0, gamma2=3.0/5.0, gamma3=3.0/10.0;
   double beta1, beta2, beta3;
   double u1, u2, u3;
   double w1, w2, w3;
   
   beta1 = (13.0/12.0)*pow((um2 - 2.0*um1 + u0),2) +
   (1.0/4.0)*pow((um2 - 4.0*um1 + 3.0*u0),2);
   beta2 = (13.0/12.0)*pow((um1 - 2.0*u0 + up1),2) +
   (1.0/4.0)*pow((um1 - up1),2);
   beta3 = (13.0/12.0)*pow((u0 - 2.0*up1 + up2),2) +
   (1.0/4.0)*pow((3.0*u0 - 4.0*up1 + up2),2);
   
   w1 = gamma1 / pow(eps+beta1, 2);
   w2 = gamma2 / pow(eps+beta2, 2);
   w3 = gamma3 / pow(eps+beta3, 2);
   
   u1 = (1.0/3.0)*um2 - (7.0/6.0)*um1 + (11.0/6.0)*u0;
   u2 = -(1.0/6.0)*um1 + (5.0/6.0)*u0 + (1.0/3.0)*up1;
   u3 = (1.0/3.0)*u0 + (5.0/6.0)*up1 - (1.0/6.0)*up2;
   
   return (w1 * u1 + w2 * u2 + w3 * u3)/(w1 + w2 + w3);
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
double shock_ind(double a, double b, double c)
{
   return fabs(a - 2.0*b + c) / fabs(a + 2.0*b + c);
   //static const double om = 1.0;
   //return fabs(a - 2.0*b + c) / (om*(fabs(a-b) + fabs(b-c)) + (1-om)*fabs(a + 2.0*b + c));
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
      void reconstruct (const unsigned int face,
                        vector<double>& left,
                        vector<double>& right) const;               
      void ent_diss_flux(const vector<double>& left,
                                 const vector<double>& right,
                                 double rho,
                                 double u,
                                 double a,
                                 double betal,
                                 double betar,
                                 vector<double>& flux) const;
      void ent_diss_flux_2(const vector<double>& left_m1,
                              const vector<double>& left,
                              const vector<double>& right,
                              const vector<double>& right_p1,
                              const vector<double>& right_p2,
                              double rho_m1,
                              double u_m1,
                              double a_m1,
                              double betal_m1,
                              double betar_m1,
                              double rho,
                              double u,
                              double a,
                              double betal,
                              double betar,
                              double rho_p1,
                              double u_p1,
                              double a_p1,
                              double betal_p1,
                              double betar_p1,
                              vector<double>& flux) const  ;                       
      void keps_ent_flux(const vector<double>& left,
                         const vector<double>& right,
                         vector<double>& flux) const;
      void keps2_ent_flux(const int& f,
                          const vector<double>& left,
                          const vector<double>& right,
                          vector<double>& flux) const;
      void roe_flux(const vector<double>& left,
                         const vector<double>& right,
                         vector<double>& flux) const;
      void roe_fixed_flux(const vector<double>& left,
                         const vector<double>& right,
                         vector<double>& flux) const;                   
      void rusanov_flux ( const vector<double>& left,
                         const vector<double>& right,
                         vector<double>& flux) const;                             
      void num_flux (const int&  f,
                     const vector<double>&,
                     const vector<double>&,
                     vector<double>&) const;              
      void compute_face_values ();
      void compute_residual ();
      void compute_residual_norm ();
      void update_solution (const unsigned int rk);
      void update_ent_res(const unsigned int rk);
      double entropy( const vector<double>& prim_var) const;                                            
      void output ();
      double conval(int i, int j) const;
      
      FluxScheme flux_scheme;
      double kepes_diss;
      ReconstructionScheme recon_scheme;
      TimeIntegrationScheme time_scheme;

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
      int    max_iter;
      vector<double> xc;
      vector<double> xf;
   
      vector< vector<double> > primitive;
      vector< vector<double> > primitive_old;
      vector< vector<double> > residual;
      vector< vector<double> > conserved;
      vector< vector<double> > conserved_old;
      vector< vector<double> > conl, conr;
      vector<double> res_norm;
      
      double M; // Defining Mach numbers for the initial conditions
                 // of certain test cases
      int itermod; 
      int test_case;
      
      int ent_diss_order; // order of entropy based matrix dissiaption
};

//------------------------------------------------------------------------------
// Constructor:
// INSTRUCTIONS FOR PARAMETER SELECTION
//
// Choosing time integration scheme (time_scheme):
//    RK1, SSPRK3, JAMESON_RK4
//
// Choosing flux (flux_scheme):
//    KEPSSENT     :: Exactly entropy conservative and kinetic energy preserving scheme with
//                    entropy variable based matrix dissipation
//                    We use kepes_rusanov_roe hybrid dissipation if kep_rus_roe_hyb = 1  
//                    FOR TECNO TYPE SCHEME SET ent_diss_order = 1
//    ROE          :: Original Roe scheme
//    ROEFIXED     :: Roe scheme with entropy fix
//    RUSSANOV     :: Rusanov flux 
//
//
// Scalar dissipation when flux_scheme = KEPS
//    1.0 = full dissipation, 0.0 = only density dissipation
//    K2, K4 to be specified with the test cases. These control the second and fourth order
//    dissipation respectively
//
// Choose type of reconstruction (recon_scheme):
//    FIRST,MINMOD,VANLEER,WENO
//    If using TECNO scheme, set recon_scheme to FIRST
//
// Order of matrix dissipation when used with KEPSSENT (ent_diss_order)
//    1 = first order
//    2 = second order TECNO 
//
// Eigenvalue augmentation for entropy consistency
//    alpha : for acoustic eigenvalue jump augmentation
//    beta_upwind : for acoustic upwind eigenvalue modification
//
// Set test case by choosing from one of the available test cases
//
//------------------------------------------------------------------------------
FVProblem::FVProblem ()
{   
   n_var  = 3;
   
   time_scheme = SSPRK3;
   
   if(time_scheme == RK1)
      nrk = 1;
   else if(time_scheme == SSPRK3)
      nrk = 3;
   else if(time_scheme == JAMESON_RK4)
      nrk = 4;
   else
   {
      cout<<"Unknown time integration scheme"<<endl;
      cout<<"Possible options: RK1, RK3, JAMESON_RK4"<<endl; 
   }           
   
   flux_scheme = KEPSSENT; 

   kepes_diss = 0.0;
   
   ent_diss_order = 1; 
   
   recon_scheme = FIRST;
   
   if(flux_scheme == KEPSSENT && ent_diss_order == 2 && recon_scheme != FIRST)
   {
      cout<<"Can only use FIRST order reconstruction option with TECNO scheme"<<endl;
      exit(0);
   }          
   
   alpha = 0.0/1.0; 
   beta_upwind = 0.0/6.0; 

   test_case = 11;
   max_iter = 10000000;

   if(test_case == 1)
   {
      // Sod shock tube case
      GAMMA = 1.4;
      gas_const = 1.0;
      final_time = 0.2;
      itermod = 10;
      n_cell = 100;
      cfl    = 0.2;
      K2 = 0.0;
      K4 = 0.0;

      xmin   = 0.0;
      xmax   = 1.0;
      xmid   = 0.5 * (xmin + xmax);

      d_left  = 1.0;
      d_right = 0.125;
   
      u_left  = 0.0;
      u_right = 0.0;
   
      p_left  = 1.0;
      p_right = 0.1;

   }
   else if(test_case == 11)
   {
      // modified sod
      // Entropy violation occurs for roe scheme
      GAMMA = 1.4;
      gas_const = 1.0;
      final_time = 0.2;
      itermod = 10;
      n_cell = 100;
      cfl    = 0.3;
      K2 = 0.0;
      K4 = 0.0;
      
      xmin   = 0.0;
      xmax   = 1.0;
      xmid   = 0.3;
      
      d_left  = 1.0;
      d_right = 0.125;
      
      u_left  = 0.75;
      u_right = 0.0;
      
      p_left  = 1.0;
      p_right = 0.1;
      
   }
   else if(test_case == 12)
   {
      // stationary contact
      GAMMA = 1.4;
      gas_const = 1.0;
      final_time = 100.0;
      itermod = 100;
      n_cell = 26;
      cfl    = 0.4;
      K2 = 0.0;
      K4 = 0.0;
      
      xmin   = 0.0;
      xmax   = 1.0;
      xmid   = 0.5;
      
      d_left  = 10.0;
      d_right = 1.0;
      
      u_left  = 0.0;
      u_right = 0.0;
      
      p_left  = 1.0;
      p_right = 1.0;

   }
   else if(test_case == 2)
   {
      // shock structure case
      GAMMA = 5.0/3.0;
      gas_const = 0.5;
      final_time = 200.0;
      itermod = 1000;
      n_cell = 100;
      cfl    = 0.1;
      K2 = 0.0;
      K4 = 1.0/100.0;

      xmin   = -0.25;
      xmax   =  0.00;
      xmid   = 0.5 * (xmin + xmax);
      M = 1.5;
      double M2 = pow(M, 2);

      d_left = 1.0;
      u_left = 1.0;
      p_left = 1.0/GAMMA/M2;

      d_right= (GAMMA+1.0)*M2/(2.0+(GAMMA-1.0)*M2)*d_left;
      u_right= ((GAMMA-1.0)/(GAMMA+1.0)+2.0/(GAMMA+1.0)/M2)*u_left;
      p_right= (2.0*GAMMA/(GAMMA+1.0)*M2-(GAMMA-1)/(GAMMA+1.0))*p_left;

   }
   else if(test_case == 3)
   {
      // Stationary shock
      GAMMA = 1.4;
      gas_const = 1.0;
      final_time = 5.0;
      itermod = 100;
      n_cell = 200;
      cfl    = 0.4;
      K2 = 0.0;
      K4 = 0.0/300.0;
      
      xmin   = 0.0;
      xmax   = 1.0;
      xmid   = 0.5 * (xmin + xmax);
      
      M = 20;
      
      d_left  = 1.0;
      d_right = 1.0/(2.0/(GAMMA+1)/pow(M,2) + (GAMMA-1.0)/(GAMMA+1.0));		
      
      u_left  = 1.0;
      u_right = 1.0/d_right;
      
      p_left  = 1.0/(GAMMA*pow(M,2));
      p_right = p_left * (2.0*GAMMA*pow(M,2)/(GAMMA+1.0) - (GAMMA-1.0)/(GAMMA+1.0));
      
   }
   else if(test_case == 4)
   {
      // Rarefaction
      GAMMA = 1.4;
      gas_const = 1.0;
      final_time = 0.1;
      itermod = 1;
      n_cell = 50;
      cfl    = 0.2;
      K2 = 0.0;
      K4 = 0.0;
      
      xmin   = 0.0;
      xmax   = 1.0;
      xmid   = 0.5 * (xmin + xmax);
      
      M = 3.0;
      
      d_right  = 1.0;
      d_left = 1.0/(2.0/(GAMMA+1)/pow(M,2) + (GAMMA-1.0)/(GAMMA+1.0));
      
      u_right  = 1.0;
      u_left = 1.0/d_left;
      
      p_right  = 1.0/(GAMMA*pow(M,2));
      p_left = p_right * (2.0*GAMMA*pow(M,2)/(GAMMA+1.0) - (GAMMA-1.0)/(GAMMA+1.0));
      
   }
   else if(test_case == 101)
   {
      // Lax Problem
      GAMMA = 5.0/3.0;
      gas_const = 1.0;
      final_time = 0.3;
      itermod = 10;
      n_cell = 100;
      cfl    = 0.4;
      K2 = 0.0;
      K4 = 0.0;
      
      xmin   = 0.0;
      xmax   = 3.0;
      xmid   = 1.5;
      
      d_left  = 0.445;
      d_right = 0.5;
      
      u_left  = 0.698;
      u_right = 0.0;
      
      p_left  = 3.528;
      p_right = 0.571;
      
   }
   else if(test_case == 102)
   {
      // Low Density Problem
      GAMMA = 1.4;
      gas_const = 1.0;
      final_time = 0.12;
      itermod = 100;
      n_cell =5000;
      cfl    = 0.6;
      K2 = 0.0;
      K4 = 0.0;
      
      xmin   = 0.0;
      xmax   = 1.0;
      xmid   = 0.5;
      
      d_left  = 1.0;
      d_right = 1.0;
      
      u_left  = -2.0;
      u_right = 2.0;
      
      p_left  = 0.4;
      p_right = 0.4;
      
   }
   else if(test_case == 103)
   {
      // Sod & Lax Problem
      GAMMA = 1.4;
      gas_const = 1.0;
      final_time = 0.09;
      itermod = 5;
      n_cell = 50;
      cfl    = 0.4;
      K2 = 0.0;
      K4 = 0.0;
      
      xmin   = 0.0;
      xmax   = 1.0;
      xmid   = 0.5;
      
      d_left  = 3.857;
      d_right = 1.0;
      
      u_left  = 0.92;
      u_right = 3.55;
      
      p_left  = 10.333;
      p_right = 1.0;
      
   }
   else if(test_case == 104)
   {
      // Slow moving weak shock Problem
      GAMMA = 1.4;
      gas_const = 1.0;
      final_time = 1.0;
      itermod = 50;
      n_cell = 100;
      cfl    = 0.1;
      K2 = 0.0;
      K4 = 0.0;
      
      xmin   = 0.0;
      xmax   = 1.0;
      xmid   = 0.5;
      
      //d_left  = 1;
      //d_right = 0.9275;
      
      //u_left  = -1.0;
      //u_right = -1.0781;
      
      //p_left  = 1.0;
      //p_right = 0.9;
      
      d_left  = 3.86;
      d_right = 1.0;
      
      u_left  = -0.81;
      u_right = -3.44;
      
      p_left  = 10.3300108;
      p_right = 1.0;
      
   }
   
     else if(test_case == 105)
   {
      //Noh's  Problem
      GAMMA = 5.0/3.0;
      gas_const = 1.0;
      final_time = 1.0;
      itermod = 10;
      n_cell = 200;
      cfl    = 0.1;
      K2 = 0.0;
      K4 = 0.0;
      
      M = 100;
      double M2 = pow(M,2);      

      xmin   = 0.0;
      xmax   = 1.0;
      xmid   = 0.5;
      
      d_left  = 1.0;
      d_right = 1.0;
      
      u_left  = 1.0;
      u_right = -1.0;
      
      p_left  = 1.0/GAMMA/M2;
      p_right = p_left;
      
   } 

   else if(test_case == 106)
   {
      //Toro Test 3
      GAMMA = 1.4;
      gas_const = 1.0;
      final_time = 0.012;
      itermod = 100;
      n_cell = 140;
      cfl    = 0.1;
      K2 = 0.0;
      K4 = 0.0;     

      xmin   = 0.0;
      xmax   = 1.4;
      xmid   = 0.7;
      
      d_left  = 1.0;
      d_right = 1.0;
      
      u_left  = 0.0;
      u_right = 0.0;
      
      p_left  = 1000.0;
      p_right = 0.01;
      
   } 
   
   else if(test_case == 107)
   {
      //Toro Test 4
      GAMMA = 1.4;
      gas_const = 1.0;
      final_time = 0.035;
      itermod = 100;
      n_cell = 400;
      cfl    = 0.2;
      K2 = 0.0;
      K4 = 0.0;     

      xmin   = 0.0;
      xmax   = 1.0;
      xmid   = 0.4;
      
      d_left  = 5.99924;
      d_right = 5.99242;
      
      u_left  = 19.5975;
      u_right = -6.19633;
      
      p_left  = 460.894;
      p_right = 46.0950;
      
   }
   
   else if(test_case == 108)
   {
      //Toro Test 5: Slow moving contact
      GAMMA = 1.4;
      gas_const = 1.0;
      final_time = 0.01;
      itermod = 100;
      n_cell = 800;
      cfl    = 0.1;
      K2 = 0.0;
      K4 = 0.0;     

      xmin   = 0.0;
      xmax   = 2.0;
      xmid   = 1.0;
      
      d_left  = 1.0;
      d_right = 1.0;
      
      u_left  = -19.59745;
      u_right = -19.59745;
      
      p_left  = 1000.0;
      p_right = 0.01;
      
   }
   
   else if(test_case == 109)
   {
      // Toro Test 7: moving isolated contact
      GAMMA = 1.4;
      gas_const = 1.0;
      final_time = 2.0;
      itermod = 10;
      n_cell = 26;
      cfl    = 0.4;
      K2 = 0.0;
      K4 = 0.0;
      
      xmin   = 0.0;
      xmax   = 1.0;
      xmid   = 0.5;
      
      d_left  = 1.4;
      d_right = 1.0;
      
      u_left  = 0.1;
      u_right = 0.1;
      
      p_left  = 1.0;
      p_right = 1.0;
      
   }
   
   else if(test_case == 111)
   {
      // Lax Problem (S. MISHRA'S IC)
      GAMMA = 1.4;
      gas_const = 1.0;
      final_time = 1.3;
      itermod = 10;
      n_cell = 100;
      cfl    = 0.4;
      K2 = 0.0;
      K4 = 0.0;
      
      xmin   = 0.0;
      xmax   = 10.0;
      xmid   = 5.0;
      
      d_left  = 0.445;
      d_right = 0.5;
      
      u_left  = 0.698;
      u_right = 0.0;
      
      p_left  = 3.528;
      p_right = 0.571;
      
   }
   
   else if(test_case == 112)
   {
      // Stationary soln (Fredrich Coquel)
      GAMMA = 1.4;
      gas_const = 1.0;
      final_time = 0.2;
      itermod = 10;
      n_cell = 50;
      cfl    = 0.4;
      K2 = 0.0;
      K4 = 0.0;
      
      xmin   = 0.0;
      xmax   = 1.0;
      xmid   = 0.5;
      
      
      double T_left = 2442;
      double T_right = 346;
      
      p_left  = 650000;
      p_right = 1000;
      
      d_left  = p_left/T_left/gas_const;
      d_right = p_right/T_right/gas_const;
      
      u_left  = 0.0;
      u_right = 0.0;
      
   }
   else
   {
      cout<<"Unknown Test Case"<<endl;
      exit(0);
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
   primitive_old.resize (n_cell, vector<double>(n_var));
   residual.resize (n_cell, vector<double>(n_var));
   conserved.resize (n_cell, vector<double>(n_var));
   conserved_old.resize (n_cell, vector<double>(n_var));
   conl.resize (n_cell, vector<double>(n_var));
   conr.resize (n_cell, vector<double>(n_var));

   res_norm.resize (n_var);
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
// return conserved variable value
//------------------------------------------------------------------------------
double FVProblem::conval (int i, int j) const
{
   if(i>=0 && i<=n_cell-1) 
      return conserved[i][j];
   else if(i<0) 
      return conserved[0][j];
   else 
      return conserved[n_cell-1][j];
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
   dt *= cfl;
   
}

//------------------------------------------------------------------------------
// Convert conserved to primitive
//------------------------------------------------------------------------------
void con2prim (vector<double>& con, vector<double>& prim)
{
   prim[0] = con[0];
   prim[1] = con[1]/con[0];
   prim[2] = (GAMMA-1.0) * (con[2] - 0.5 * pow(con[1], 2.0) / con[0]);
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
// Reconstruct left/right state at a face
//------------------------------------------------------------------------------
void FVProblem::reconstruct (const unsigned int face,
                             vector<double>& left,
                             vector<double>& right) const
{
   if(recon_scheme == FIRST)
   {
      left = primitive[face-1];
      right= primitive[face];
   }
   else if(recon_scheme == MINMOD || recon_scheme == VANLEER)
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
   else if(recon_scheme == WENO)
   {
      vector<double> cl(n_var), cr(n_var);
      for(unsigned int j=0; j<n_var; ++j)
      {
         cl[j] = weno5(conval(face-3,j), 
                       conval(face-2,j), 
                       conval(face-1,j), 
                       conval(face,j), 
                       conval(face+1,j));
         cr[j] = weno5(conval(face+2,j), 
                       conval(face+1,j), 
                       conval(face,j), 
                       conval(face-1,j), 
                       conval(face-2,j));
      }
      con2prim(cl, left);
      con2prim(cr, right);
   }
   else
   {
      cout<<"Unknown type of variable reconstruction specified."<<endl; 
      cout<<"Possible options FIRST, MINMOD, VANLEER, WENO"<<endl;
      exit(0);
   }
}


//------------------------------------------------------------------------------
// Numerical flux function
// exactly entropy consistent
//------------------------------------------------------------------------------
void FVProblem::keps2_ent_flux(const int& f,
                                const vector<double>& left,
                               const vector<double>& right,
                               vector<double>& flux) const
{   
   vector<double> left_m1(n_var), right_p1(n_var), right_p2(n_var);
   if(f==0)
   {
      left_m1  = prim_left;
      right_p1 = primitive[f+1];
      right_p2 = primitive[f+2];
   }
   else if(f==1)
   {
      left_m1  = prim_left;
      right_p1 = primitive[f+1];
      right_p2 = primitive[f+2];
   }
   else if(f==n_face-3)
   {
      left_m1  = primitive[f-2];
      right_p1 = primitive[f+1];
      right_p2 = prim_right;
   }
   else if(f==n_face-2)
   {
      left_m1  = primitive[f-2];
      right_p1 = prim_right;
      right_p2 = prim_right;
   }
   else if(f==n_face-1)
   {
      left_m1  = primitive[f-2];
      right_p1 = prim_right;
      right_p2 = prim_right;
   }
   else
   {
      left_m1  = primitive[f-2];
      right_p1 = primitive[f+1];
      right_p2 = primitive[f+2];
   }
 
   double rho = logavg (left[0], right[0]);
   double u   = 0.5 * (left[1] + right[1]);
   double u2  = 0.5 * (pow(left[1],2) + pow(right[1],2));
   double Tl  = temperature (left);
   double Tr  = temperature (right);
   double betal = 1.0/(2.0*gas_const*Tl);
   double betar = 1.0/(2.0*gas_const*Tr);
   double beta  = logavg(betal, betar);
   double rho_a = 0.5*(left[0]+right[0]);
   double beta_a = 0.5*(betal+betar);
   double p   = 0.5 * rho_a / beta_a;
    
   flux[0] = rho * u;
   flux[1] = p + u * flux[0];
   flux[2] = (1.0/(2.0*(GAMMA-1.0)*beta) - 0.5*u2) * flux[0] + u * flux[1];
       
   double a = sqrt(0.5 * GAMMA / beta);   
       
     
   // Add entropy dissipation
   if(ent_diss_order == 1)
      ent_diss_flux(left, right, rho, u, a, betal, betar, flux);
   else if(ent_diss_order == 2)
   {
	  double rho_m1 = logavg (left_m1[0], left[0]);
      double u_m1   = 0.5 * (left_m1[1] + left[1]);
      double Tl_m1  = temperature (left_m1);
      double Tr_m1  = temperature (left);
      double betal_m1 = 1.0/(2.0*gas_const*Tl_m1);
      double betar_m1 = 1.0/(2.0*gas_const*Tr_m1);
      double beta_m1  = logavg(betal_m1, betar_m1);
       
      double a_m1 = sqrt(0.5 * GAMMA / beta_m1); 
   
      double rho_p1 = logavg (right_p1[0], right_p2[0]);
      double u_p1   = 0.5 * (right_p1[1] + right_p2[1]);
      //double u2_p1  = 0.5 * (pow(right_p1[1],2) + pow(right_p2[1],2));
      double Tl_p1  = temperature (right_p1);
      double Tr_p1  = temperature (right_p2);
      double betal_p1 = 1.0/(2.0*gas_const*Tl_p1);
      double betar_p1 = 1.0/(2.0*gas_const*Tr_p1);
      double beta_p1  = logavg(betal_p1, betar_p1);
       
      double a_p1 = sqrt(0.5 * GAMMA / beta_p1); 
     
	  ent_diss_flux_2(left_m1,left, right, right_p1,right_p2, 
	                rho_m1, u_m1, a_m1, betal_m1, betar_m1, 
	                rho, u, a, betal, betar, 
	                rho_p1, u_p1,  a_p1, betal_p1, betar_p1,  
	                flux); 
   }
   
}

//------------------------------------------------------------------------------
// Entropy variable matrix dissipation flux
// We use kepes_rusanov_roe hybrid dissipation if kep_rus_roe_hyb = 1
//------------------------------------------------------------------------------
void FVProblem::ent_diss_flux(const vector<double>& left,
                              const vector<double>& right,
                              double rho,
                              double u,
                              double a,
                              double betal,
                              double betar,
                              vector<double>& flux) const
{
   // Add entropy dissipation
   //double a = sqrt(0.5 * GAMMA / beta);
   double H = a*a/(GAMMA-1.0) + 0.5*u*u;
   double R[3][3];
   R[0][0] = R[0][1] = R[0][2] = 1.0;
   R[1][0] = u-a; R[1][1] = u; R[1][2] = u + a;
   R[2][0] = H - u * a; R[2][1] = 0.5*u*u; R[2][2] = H + u * a;
   
   double ul = left[1];
   double ur = right[1];
   double al = sqrt(0.5 * GAMMA / betal);
   double ar = sqrt(0.5 * GAMMA / betar);
   
   double LambdaL[] = { ul-al, ul, ul+al };
   double LambdaR[] = { ur-ar, ur, ur+ar };
   
   double l1,l2,l3;
   
   l1 = fabs(u-a);
   l2 = fabs(u);
   l3 = fabs(u+a);
   
   double phi_switch = 0.0;
   
   double LambdaRoe[] = { (1+beta_upwind)*l1 + alpha*fabs(LambdaR[0]-LambdaL[0]), 
                           l2, 
                           (1+beta_upwind)*l3 + alpha*fabs(LambdaR[2]-LambdaL[2])
                         };
                     
   double LambdaRus[] = { fabs(u) + a, fabs(u) + a, fabs(u) + a};
   
   double Lambda[] = { (1.0 - phi_switch)*LambdaRoe[0] + phi_switch*LambdaRus[0],
                        (1.0 - phi_switch)*LambdaRoe[1] + phi_switch*LambdaRus[1],
                        (1.0 - phi_switch)*LambdaRoe[2] + phi_switch*LambdaRus[2]
                      };  
   
   double S[] = { 0.5*rho/GAMMA, (GAMMA-1.0)*rho/GAMMA, 0.5*rho/GAMMA };
   double D[] = { Lambda[0]*S[0], Lambda[1]*S[1], Lambda[2]*S[2] };
   
   double sl = log(left[2]) - GAMMA * log(left[0]);
   double sr = log(right[2]) - GAMMA * log(right[0]);
   double ds = sr - sl;
   
   // Jump in entropy variables
   double dv[] = { - ds/(GAMMA-1.0) - (betar*right[1]*right[1] - betal*left[1]*left[1]),
                     2.0*(betar*right[1] - betal*left[1]), 
                    -2.0*(betar-betal) };
   
   double Diff[] = {0, 0, 0};
   for(unsigned int i=0; i<3; ++i)
   {
      for(unsigned int j=i; j<3; ++j)
      {
		 double RDRT = 0;
		 for(unsigned int k=0; k<3; ++k) 
		    RDRT += R[i][k]*D[k]*R[j][k];
		 Diff[i] += RDRT*dv[j];
		 if(i!=j)
		    Diff[j] += RDRT*dv[i];
	  }
	  flux[i] -= 0.5*Diff[i];
   }
   
}


//------------------------------------------------------------------------------
// Entropy variable matrix dissipation flux 2nd order
// We use kepes_rusanov_roe hybrid dissipation if kep_rus_roe_hyb = 1
//------------------------------------------------------------------------------
void FVProblem::ent_diss_flux_2(const vector<double>& left_m1,
                              const vector<double>& left,
                              const vector<double>& right,
                              const vector<double>& right_p1,
                              const vector<double>& right_p2,
                              double rho_m1,
                              double u_m1,
                              double a_m1,
                              double betal_m1,
                              double betar_m1,
                              double rho,
                              double u,
                              double a,
                              double betal,
                              double betar,
                              double rho_p1,
                              double u_p1,
                              double a_p1,
                              double betal_p1,
                              double betar_p1,
                              vector<double>& flux) const
{
   // Add entropy dissipation
   //double a_m1 = sqrt(0.5 * GAMMA / beta_m1);
   double H_m1 = a_m1*a_m1/(GAMMA-1.0) + 0.5*u_m1*u_m1;
   double R_m1[3][3];
   R_m1[0][0] = R_m1[0][1] = R_m1[0][2] = 1.0;
   R_m1[1][0] = u_m1-a_m1; R_m1[1][1] = u_m1; R_m1[1][2] = u_m1 + a_m1;
   R_m1[2][0] = H_m1 - u_m1 * a_m1; R_m1[2][1] = 0.5*u_m1*u_m1; R_m1[2][2] = H_m1 + u_m1 * a_m1;
   
   //double a = sqrt(0.5 * GAMMA / beta);
   double H = a*a/(GAMMA-1.0) + 0.5*u*u;
   double R[3][3];
   R[0][0] = R[0][1] = R[0][2] = 1.0;
   R[1][0] = u-a; R[1][1] = u; R[1][2] = u + a;
   R[2][0] = H - u * a; R[2][1] = 0.5*u*u; R[2][2] = H + u * a;
   
   //double a_p1 = sqrt(0.5 * GAMMA / beta_p1);
   double H_p1 = a_p1*a_p1/(GAMMA-1.0) + 0.5*u_p1*u_p1;
   double R_p1[3][3];
   R_p1[0][0] = R_p1[0][1] = R_p1[0][2] = 1.0;
   R_p1[1][0] = u_p1-a_p1; R_p1[1][1] = u_p1; R_p1[1][2] = u_p1 + a_p1;
   R_p1[2][0] = H_p1 - u_p1 * a_p1; R_p1[2][1] = 0.5*u_p1*u_p1; R_p1[2][2] = H_p1 + u_p1 * a_p1;
   
   double ul = left[1];
   double ur = right[1];
   double al = sqrt(0.5 * GAMMA / betal);
   double ar = sqrt(0.5 * GAMMA / betar);
   
   double LambdaL[] = { ul-al, ul, ul+al };
   double LambdaR[] = { ur-ar, ur, ur+ar };
   
   double l1,l2,l3;
   
   l1 = fabs(u-a);
   l2 = fabs(u);
   l3 = fabs(u+a);
   //l1 = l3 = min(l1, l3);
   //l1 = l3 = max(l1, l3);
   //l1 = l3 = 0.5*(l1 + l3);
   //l1 = l3 = sqrt(l1*l3);
   //l1 = l3 = 2*l1*l3/(l1+l3);
   
   
   double phi_switch = 0.0;
   
   double LambdaRoe[] = { (1+beta_upwind)*l1 + alpha*fabs(LambdaR[0]-LambdaL[0]), 
                           l2, 
                           (1+beta_upwind)*l3 + alpha*fabs(LambdaR[2]-LambdaL[2])
                         };
                     
   double LambdaRus[] = { fabs(u) + a, fabs(u) + a, fabs(u) + a};
   
   double Lambda[] = { (1.0 - phi_switch)*LambdaRoe[0] + phi_switch*LambdaRus[0],
                        (1.0 - phi_switch)*LambdaRoe[1] + phi_switch*LambdaRus[1],
                        (1.0 - phi_switch)*LambdaRoe[2] + phi_switch*LambdaRus[2]
                      };  
   
   double S[] = { 0.5*rho/GAMMA, (GAMMA-1.0)*rho/GAMMA, 0.5*rho/GAMMA };
   double D[] = { Lambda[0]*S[0], Lambda[1]*S[1], Lambda[2]*S[2] };
   
   
   double sl_m1 = log(left_m1[2]) - GAMMA * log(left_m1[0]);
   double sr_m1 = log(left[2]) - GAMMA * log(left[0]);
   double ds_m1 = sr_m1 - sl_m1;
   
   double sl = log(left[2]) - GAMMA * log(left[0]);
   double sr = log(right[2]) - GAMMA * log(right[0]);
   double ds = sr - sl;
   
   double sl_p1 = log(right[2]) - GAMMA * log(right[0]);
   double sr_p1 = log(right_p1[2]) - GAMMA * log(right_p1[0]);
   double ds_p1 = sr_p1 - sl_p1;
   
   // Jump in entropy variables
   double dv_m1[] = { - ds_m1/(GAMMA-1.0) - (betar_m1*left[1]*left[1] - betal_m1*left_m1[1]*left_m1[1]),
                     2.0*(betar_m1*left[1] - betal_m1*left_m1[1]), 
                    -2.0*(betar_m1-betal_m1) };
                    
   double dv[] = { - ds/(GAMMA-1.0) - (betar*right[1]*right[1] - betal*left[1]*left[1]),
                     2.0*(betar*right[1] - betal*left[1]), 
                    -2.0*(betar-betal) };
                    
   double dv_p1[] = { - ds_p1/(GAMMA-1.0) - (betar_p1*right_p1[1]*right_p1[1] - betal_p1*right[1]*right[1]),
                     2.0*(betar_p1*right_p1[1] - betal_p1*right[1]), 
                    -2.0*(betar_p1-betal_p1) };                 
   
   double dw_m1[3] , dw[3], dw_p1[3];
   
   for(unsigned int i=0; i<3; ++i)
   {
	  dw_m1[i] = 0.0; 
      dw[i] = 0.0;
	  dw_p1[i] = 0.0; 
      for(unsigned int j=0; j<3; ++j)
      {
		  dw_m1[i] += R_m1[j][i] * dv_m1[j]; 
		  dw[i] += R[j][i] * dv[j];
		  dw_p1[i] += R_p1[j][i] * dv_p1[j];
	  }
   }	  
	  	  
   double dw_tilda[3];
   
   for(unsigned int i=0; i<3; ++i)
   {
		dw_tilda[i]= (1 - 0.5*(minmod2(dw_m1[i],dw[i]) + minmod2(dw_p1[i],dw[i])))*dw[i]; 
   }
   		  
   
   // LSW = L*S*dWtilda
   double LSW[3];
   for(unsigned int i=0; i<3; ++i)
      LSW[i] = D[i]*dw_tilda[i];
   
   for(unsigned int i=0; i<3; ++i)
   {
      double Diff = 0.0;
      for(unsigned int j=0; j<3; ++j)
          Diff += R[i][j]*LSW[j];
      
      flux[i] -= 0.5*Diff;
   }
   
}

//------------------------------------------------------------------------------
// Numerical flux function
// Original Roe Scheme
//------------------------------------------------------------------------------
void FVProblem::roe_flux(const vector<double>& left,
                         const vector<double>& right,
                         vector<double>& flux) const
{   
   double fl  = sqrt(left[0]);
   double fr  = sqrt(right[0]);
   double u   = (fl*left[1] + fr*right[1])/(fl + fr);
   double Tl  = temperature (left);
   double Tr  = temperature (right);

   double Hl  = GAMMA*gas_const*Tl/(GAMMA-1.0) + 0.5*pow(left[1],2);
   double Hr  = GAMMA*gas_const*Tr/(GAMMA-1.0) + 0.5*pow(right[1],2);
   
   double El  = left[2]/(GAMMA-1.0) + 0.5*left[0]*pow(left[1],2);
   double Er  = right[2]/(GAMMA-1.0) + 0.5*right[0]*pow(right[1],2);
   
   // average of fluxes
   flux[0] = 0.5*(left[0]*left[1] + right[0]*right[1]);
   flux[1] = 0.5*(left[2] + left[0] * pow(left[1],2) + 
                  right[2] + right[0] * pow(right[1],2));
   flux[2] = 0.5*(Hl*left[0]*left[1] + Hr*right[0]*right[1]);
   
   
   // Add conservative dissipation
   double H = (fl*Hl + fr*Hr)/(fl + fr);
   double a = sqrt((GAMMA-1.0)*(H - 0.5*u*u));
   double R[3][3];
   R[0][0] = R[0][1] = R[0][2] = 1.0;
   R[1][0] = u-a; R[1][1] = u; R[1][2] = u + a;
   R[2][0] = H - u * a; R[2][1] = 0.5*u*u; R[2][2] = H + u * a;
   
   double ul = left[1];
   double ur = right[1];
   
   double Lambda[] = { fabs(u-a), fabs(u), fabs(u+a)};
   
   double dU[] = {
      right[0] - left[0], right[0]*ur - left[0]*ul, Er - El
   };
   
   double aa[3];
   aa[1] = (GAMMA-1.0)/(a*a) * (dU[0]*(H-u*u) + u*dU[1] - dU[2]);
   aa[0] = 0.5/a * (dU[0]*(u+a) - dU[1] - a * aa[1]);
   aa[2] = dU[0] - aa[0] - aa[1];
   
   for(unsigned int i=0; i<3; ++i)
      for(unsigned int j=0; j<3; ++j)
         flux[i] -= 0.5 * aa[j] * Lambda[j] * R[i][j];
   
}

//------------------------------------------------------------------------------
// Roe Scheme with entropy fix
//------------------------------------------------------------------------------
void FVProblem::roe_fixed_flux(const vector<double>& left,
                         const vector<double>& right,
                         vector<double>& flux) const
{   
   double fl  = sqrt(left[0]);
   double fr  = sqrt(right[0]);
   double u   = (fl*left[1] + fr*right[1])/(fl + fr);
   double Tl  = temperature (left);
   double Tr  = temperature (right);

   double Hl  = GAMMA*gas_const*Tl/(GAMMA-1.0) + 0.5*pow(left[1],2);
   double Hr  = GAMMA*gas_const*Tr/(GAMMA-1.0) + 0.5*pow(right[1],2);
   
   double El  = left[2]/(GAMMA-1.0) + 0.5*left[0]*pow(left[1],2);
   double Er  = right[2]/(GAMMA-1.0) + 0.5*right[0]*pow(right[1],2);
   
   // average of fluxes
   flux[0] = 0.5*(left[0]*left[1] + right[0]*right[1]);
   flux[1] = 0.5*(left[2] + left[0] * pow(left[1],2) + 
                  right[2] + right[0] * pow(right[1],2));
   flux[2] = 0.5*(Hl*left[0]*left[1] + Hr*right[0]*right[1]);
   
   
   // Add conservative dissipation
   double H = (fl*Hl + fr*Hr)/(fl + fr);
   double a = sqrt((GAMMA-1.0)*(H - 0.5*u*u));
   double R[3][3];
   R[0][0] = R[0][1] = R[0][2] = 1.0;
   R[1][0] = u-a; R[1][1] = u; R[1][2] = u + a;
   R[2][0] = H - u * a; R[2][1] = 0.5*u*u; R[2][2] = H + u * a;
   
   double ul = left[1];
   double ur = right[1];
   
   // Entropy fix
   double delta = 0.1*a;
   double Lambda1, Lambda3;
   if (fabs(u-a) > delta)
		Lambda1 = fabs(u-a);
   else
        Lambda1 = (pow(u-a,2) + pow(delta,2))/(2*delta);
   if (fabs(u+a) > delta)
		Lambda3 = fabs(u+a);
   else
        Lambda3 = (pow(u+a,2) + pow(delta,2))/(2*delta);
   	
   double Lambda[] = { Lambda1, fabs(u), Lambda3};
   
   double dU[] = {
      right[0] - left[0], right[0]*ur - left[0]*ul, Er - El
   };
   
   double aa[3];
   aa[1] = (GAMMA-1.0)/(a*a) * (dU[0]*(H-u*u) + u*dU[1] - dU[2]);
   aa[0] = 0.5/a * (dU[0]*(u+a) - dU[1] - a * aa[1]);
   aa[2] = dU[0] - aa[0] - aa[1];
   
   for(unsigned int i=0; i<3; ++i)
      for(unsigned int j=0; j<3; ++j)
         flux[i] -= 0.5 * aa[j] * Lambda[j] * R[i][j];
   
}

//------------------------------------------------------------------------------
// Numerical flux function
// Rusanov Scheme
//------------------------------------------------------------------------------
void FVProblem::rusanov_flux(const vector<double>& left,
                         const vector<double>& right,
                         vector<double>& flux) const
{   
   double fl  = sqrt(left[0]);
   double fr  = sqrt(right[0]);
   double u   = (fl*left[1] + fr*right[1])/(fl + fr);
   double Tl  = temperature (left);
   double Tr  = temperature (right);

   double Hl  = GAMMA*gas_const*Tl/(GAMMA-1.0) + 0.5*pow(left[1],2);
   double Hr  = GAMMA*gas_const*Tr/(GAMMA-1.0) + 0.5*pow(right[1],2);
   
   double El  = left[2]/(GAMMA-1.0) + 0.5*left[0]*pow(left[1],2);
   double Er  = right[2]/(GAMMA-1.0) + 0.5*right[0]*pow(right[1],2);
   
   // average of fluxes
   flux[0] = 0.5*(left[0]*left[1] + right[0]*right[1]);
   flux[1] = 0.5*(left[2] + left[0] * pow(left[1],2) + 
                  right[2] + right[0] * pow(right[1],2));
   flux[2] = 0.5*(Hl*left[0]*left[1] + Hr*right[0]*right[1]);
   
   
   // Add Rusanov dissipation
   double H = (fl*Hl + fr*Hr)/(fl + fr);
   double a = sqrt((GAMMA-1.0)*(H - 0.5*u*u));
   double R[3][3];
   R[0][0] = R[0][1] = R[0][2] = 1.0;
   R[1][0] = u-a; R[1][1] = u; R[1][2] = u + a;
   R[2][0] = H - u * a; R[2][1] = 0.5*u*u; R[2][2] = H + u * a;
   
   double ul = left[1];
   double ur = right[1];
   
   double MaxLambda = fabs(u) + a;
   double Lambda[] = { MaxLambda,MaxLambda,MaxLambda};
   
   double dU[] = {
      right[0] - left[0], right[0]*ur - left[0]*ul, Er - El
   };
   
   double aa[3];
   aa[1] = (GAMMA-1.0)/(a*a) * (dU[0]*(H-u*u) + u*dU[1] - dU[2]);
   aa[0] = 0.5/a * (dU[0]*(u+a) - dU[1] - a * aa[1]);
   aa[2] = dU[0] - aa[0] - aa[1];
   
   for(unsigned int i=0; i<3; ++i)
      for(unsigned int j=0; j<3; ++j)
         flux[i] -= 0.5 * aa[j] * Lambda[j] * R[i][j];
   
}


//------------------------------------------------------------------------------
// Numerical flux function
//------------------------------------------------------------------------------
void FVProblem::num_flux(const int&            f,
                         const vector<double>& left,
                         const vector<double>& right,
                         vector<double>&       flux) const
{
   
   switch(flux_scheme)
   {
         case KEPSSENT:
            keps2_ent_flux(f, left,right,flux);
            break;
         case ROE:
            roe_flux(left, right, flux);
            break;
         case ROEFIXED:
            roe_fixed_flux(left, right, flux);
            break; 
         case RUSANOV:
            rusanov_flux(left, right, flux);
            break;       
             
   }
     
}

//------------------------------------------------------------------------------
// Compute finite volume residual
//------------------------------------------------------------------------------
void FVProblem::compute_residual ()
{
   for(unsigned int i=0; i<n_cell; ++i)
   {   for(unsigned int j=0; j<n_var; ++j)
         residual[i][j] = 0.0;
   }
   
   vector<double> flux (n_var);
   vector<double> left (n_var);
   vector<double> right(n_var); 
   
   // Flux through left boundary
   num_flux (0, prim_left, primitive[0], flux);
   for(unsigned int j=0; j<n_var; ++j)
      residual[0][j] -= flux[j];
      
   // Flux through interior faces
   for(unsigned int i=1; i<n_face-1; ++i)
   {
      reconstruct (i, left, right);
      num_flux (i, left, right, flux);
      for(unsigned int j=0; j<n_var; ++j)
      {
         residual[i-1][j] += flux[j];
         residual[i][j]   -= flux[j];           
      }
   }

   // Flux through right boundary
   num_flux (n_face-1, primitive[n_cell-1], prim_right, flux);
   for(unsigned int j=0; j<n_var; ++j)
       residual[n_cell-1][j] += flux[j];
    
}

//------------------------------------------------------------------------------
// Compute finite volume residual
//------------------------------------------------------------------------------
void FVProblem::update_solution (const unsigned int rk)
{
   if (time_scheme == RK1 || time_scheme == SSPRK3)
      for(unsigned int i=0; i<n_cell; ++i)
         for(unsigned int j=0; j<n_var; ++j)
            conserved[i][j] = arks[rk] * conserved_old[i][j] +
                              brks[rk] * (conserved[i][j] - (dt/dx) * residual[i][j]);
   else if (time_scheme == JAMESON_RK4)
      for(unsigned int i=0; i<n_cell; ++i)
         for(unsigned int j=0; j<n_var; ++j)
            conserved[i][j] = conserved_old[i][j] +
                              jameson_rks[rk] * (- (dt/dx) * residual[i][j]);   
}


//------------------------------------------------------------------------------
// Compute entropy
//------------------------------------------------------------------------------
double FVProblem::entropy(const vector<double>& prim_var) const  
{
   double s,eta;
   s = log(prim_var[2]) - GAMMA*log(prim_var[0]);
   eta = -prim_var[0]*s/(GAMMA -1);
   return eta;
}


//------------------------------------------------------------------------------
// Save solution to file
//------------------------------------------------------------------------------
void FVProblem::output ()
{
   cout<<"Saving solutions in"<<endl;
   

   //cout<<"check1"<<endl;
   string filename1 = "sol";
   string extension1 = ".dat";
   string precount;
   if     (counter <= 9)    precount = "000";
   else if(counter <= 99)   precount = "00";
   else if(counter <= 999)  precount = "0";
   else if(counter <= 9999) precount = "";
   else
   {
      cout << "Writer::output: counter is too large !!!\n";
   }
   
   stringstream ss;
   ss <<precount<<counter;
   filename1 += ss.str();
   filename1 +=extension1;
      
   cout<<filename1<<endl<<endl;   
   ofstream fo(filename1.c_str());
   for(unsigned int i=0; i<n_cell; ++i)
   {
      double sonic = sqrt(GAMMA * primitive[i][2] / primitive[i][0] );
      double mach = primitive[i][1] / sonic;
      double H = sonic*sonic/(GAMMA-1) + 0.5*pow(primitive[i][1],2);
      double s = log(primitive[i][2]) - GAMMA*log(primitive[i][0]);
      fo << xc[i] << " " 
         << primitive[i][0] << " " 
         << primitive[i][1] << " "
         << primitive[i][2] << " "
         << mach << "  " << H << "  " << s << endl;
      
   }
   fo.close ();
   counter++;
}


//------------------------------------------------------------------------------
// compute norm of residual
//------------------------------------------------------------------------------
void FVProblem::compute_residual_norm ()
{
   for(unsigned int j=0; j<n_var; ++j)
      res_norm[j] = 0;
   
   for(unsigned int i=0; i<n_cell; ++i)
      for(unsigned int j=0; j<n_var; ++j)
         res_norm[j] += pow(residual[i][j],2);
   
   for(unsigned int j=0; j<n_var; ++j)
   {
      res_norm[j] /= n_cell;
      res_norm[j] = sqrt(res_norm[j]);
   }
}

//------------------------------------------------------------------------------
// Start the computations
//------------------------------------------------------------------------------
void FVProblem::run ()
{
   make_grid_and_dofs ();
   initial_condition ();
   con_to_prim ();
   counter = 0;
   output();
   double time = 0.0;
   unsigned int iter = 0;
   while (time < final_time && iter < max_iter)
   {
      conserved_old = conserved;
      primitive_old = primitive;
      compute_dt ();
      if(time+dt > final_time) dt = final_time - time;
      for(unsigned int rk=0; rk<nrk; ++rk)
      {
         compute_residual ();
         update_solution (rk);
         con_to_prim ();
      }
      time += dt;
      ++iter;
      compute_residual_norm();
      //cout << "Iter = " << iter <<" Time = "<<time<<endl;


      if(iter % itermod == 0 || time == final_time)
      {
         cout << "Iter = " << iter << " Time = " << time << endl;
         cout << "Residual norm = ";
         cout << res_norm[0] << " " << res_norm[1] << " " << res_norm[2] << endl;
         output ();
         //exit(0);
      }
   }


}

//------------------------------------------------------------------------------
// This where it all starts
//------------------------------------------------------------------------------
int main ()
{
   FVProblem fv_problem;
   fv_problem.run ();
   
   return 0;
}
