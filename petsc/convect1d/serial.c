/* Solve linear advection equation in 1d using periodic bc,
   upwind flux, and weno scheme.
   The grid is cell-centered.
*/
#include<stdlib.h>
#include<stdio.h>
#include<math.h>

const double ark[] = {0.0, 3.0/4.0, 1.0/3.0};
const double xmin = -1.0;
const double xmax = +1.0;

//------------------------------------------------------------------------------
// Initial condition
//------------------------------------------------------------------------------
double initcond(double x)
{
   if(x >= -0.8 && x <= -0.6)
      return exp(-log(2)*pow(x+0.7,2)/0.0009);
   else if(x >= -0.4 && x <= -0.2)
      return 1.0;
   else if(x >= 0.0 && x <= 0.2)
      return 1.0 - fabs(10*(x-0.1));
   else if(x>= 0.4 && x <= 0.6)
      return sqrt(1 - 100*pow(x-0.5,2));
   else
      return 0.0;
}

//------------------------------------------------------------------------------
// Weno reconstruction of Jiang-Shu
// Return left value for face between u0, up1
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
// Save solution to file on rank = 0 process
//------------------------------------------------------------------------------
void savesol(const int nx, const double dx, double* ug)
{
   int i;
   static int count = 0;

   FILE *f;
   f = fopen("sol.dat","w");
   for(i=0; i<nx; ++i)
      fprintf(f, "%e %e\n", xmin+i*dx+0.5*dx, ug[i]);
   fclose(f);
   printf("Wrote solution into sol.dat\n");
   if(count==0)
   {
      // Initial solution is copied to sol0.dat
      system("cp sol.dat sol0.dat");
      count = 1;
   }
}

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
   int i;
   int nx=200;         // no. of cells, can change via command line
   const int sw = 3;   // stencil width
   double cfl = 0.4;

   // with ghost values
   double *u = (double*) malloc((sw+nx+sw) * sizeof(double*)) + sw;

   double dx = (xmax - xmin) / nx;
   printf("nx = %d, dx = %e\n", nx, dx);
   for(i=0; i<nx; ++i)
   {
      double x = xmin + i*dx + 0.5*dx;
      u[i] = initcond(x);
   }

   // save initial condition
   savesol(nx, dx, u);

   double res[nx], uold[nx];
   double dt = cfl * dx;
   double lam= dt/dx;
   double tfinal = 2.0, t = 0.0;

   while(t < tfinal)
   {
      if(t+dt > tfinal) // adjust dt to reach Tf
      {
         dt = tfinal - t;
         lam = dt/dx;
      }
      for(int rk=0; rk<3; ++rk) // loop for rk stages
      {
         // fill ghost values
         for (i = -sw; i < 0; ++i)
            u[i] = u[i + nx];
         for (i = nx; i < nx + sw; ++i)
            u[i] = u[i - nx];

         // First stage, store solution at time level n into uold
         if(rk==0)
            for(i=0; i<nx; ++i) uold[i] = u[i];

         for(i=0; i<nx; ++i)
            res[i] = 0.0;

         // Loop over faces and compute flux
         // dx * du_i/dt + (f_{i+1/2} - f_{i-1/2}) = 0
         // dx * du_i/dt + res_i = 0
         for(i=0; i<nx+1; ++i)
         {
            // face between i-1, i
            double uleft = weno5(u[i-3],u[i-2],u[i-1],u[i],u[i+1]);
            double flux = uleft;
            if(i==0) // first face
            {
               res[i] -= flux;
            }
            else if(i==nx) // last face
            {
               res[i-1] += flux;
            }
            else // intermediate faces
            {
               res[i]   -= flux;
               res[i-1] += flux;
            }
         }

         // Update solution
         for(i=0; i<nx; ++i)
            u[i] = ark[rk]*uold[i]
                   + (1.0-ark[rk])*(u[i] - lam * res[i]);
      }

      t += dt;
      printf("t = %f\n", t);
   }

   savesol(nx, dx, u);
}
