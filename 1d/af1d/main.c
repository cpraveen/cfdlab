// Standard active flux scheme
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define min(a,b)  ( (a) < (b) ? (a) : (b) )
#define max(a,b)  ( (a) > (b) ? (a) : (b) )

#include "hat.h"

#define NPOINTS 100

const int n = NPOINTS; // number of points
const int m = NPOINTS-1; // number of cells
const double xmin = -1.0, xmax = 1.0;
const double dx = (xmax - xmin)/(n-1);
const double Tf = 2.0;
const double cfl = 0.95;
const int tolimit = 0;

// Returns location [0,1] at which parabola takes the
// value (u0+u2)/2
double avgloc(double u0, double u1, double u2)
{
   const double eps = 1.0e-13;

   double lin = u0 - 2*u1 + u2;
   if(fabs(lin) < eps)
      return 0.5;

   double x1 = (3*u0-4*u1+u2-sqrt(5*u0*u0-16*u0*u1+16*u1*u1+6*u0*u2-16*u1*u2+5*u2*u2))/(4*(u0-2*u1+u2));
   if(x1 >= -eps && x1 - 1.0 <= eps)
   {
      x1 = max(0.0, x1);
      x1 = min(1.0, x1);
      return x1;
   }
   double x2 = (3*u0-4*u1+u2+sqrt(5*u0*u0-16*u0*u1+16*u1*u1+6*u0*u2-16*u1*u2+5*u2*u2))/(4*(u0-2*u1+u2));
   if(x2 >= -eps && x2 - 1.0 <= eps)
   {
      x2 = max(0.0, x2);
      x2 = min(1.0, x2);
      return x2;
   }

   printf("avgloc: error !!! x = %e %e", x1, x2);
   exit(0);
}

// Returns mid point value by constructing parabola
double midvalue(double ul, double ur, double ua)
{
   return 1.5 * (ua - (ul + ur)/6.0);
}

// Location of extremum location
double extloc(double u0, double u2, double ua)
{
   const double eps = 1.0e-13;
   double u1 = midvalue(u0, u2, ua);
   double lin = u0 - 2*u1 + u2;
   if(fabs(lin) < eps) // linear function
   {
      return 0.0;
   }

   return (3*u0 - 4*u1 + u2)/(4*(u0 - 2*u1 + u2));
}



// Evaluate parabola at xi in [0,1]
double sol(double xi, double ul, double ur, double ua)
{
   double umid = midvalue(ul, ur, ua);
   double val = 2*(0.5-xi)*(1-xi)*ul + 4*xi*(1-xi)*umid + 2*xi*(xi-0.5)*ur;
   return val;
}

double charint(int i, double x, double *up, double *ua, double dt)
{
   int c, vl, vr;
   double xl, xr;

   // Get cell number
   if(i==0)
   {
      c = m-1;
      vl = n-2;
      vr = n-1;
      xl = xmax - dx;
      xr = xmax;
   }
   else
   {
      c = i-1;
      vl = i-1;
      vr = i;
      xl = x - dx;
      xr = x;
   }

   double x0 = xr - dt;
   double xi = (x0 - xl) / dx;

   if(tolimit == 0)
   {
      double val = sol(xi, up[vl], up[vr], ua[c]);
      //val = max(0, val);
      //val = min(1, val);
      return val;
   }

   // We need to limit perhaps.

   // Determine if there is an interior extremum
   double xt = extloc(up[vl], up[vr], ua[c]);

   // No interior extremum
   const double eps = 1.0e-13;
   if(xt < -eps || xt-1 > eps)
   {
      double val = sol(xi, up[vl], up[vr], ua[c]);
      return val;
   }
}

void savesol(double *x, double *up, double *ua)
{
   int i, j, ns=10;
   FILE *fp;

   fp = fopen("sol.dat", "w");
   for(i=0; i<m; ++i)
   {
      double xl = x[i];
      double xr = x[i+1];
      double d  = (xr - xl)/(ns-1);
      for(j=0; j<ns; ++j)
      {
         double xx = xl + j*d;
         double xi = (xx - x[i])/dx;
         fprintf(fp,"%e %e\n", xx, sol(xi,up[i],up[i+1],ua[i]));
      }
      fprintf(fp,"\n");
   }
   fclose(fp);
}

int main()
{

   double x[n], up0[n], up1[n], up2[n], ua[m];

   double dt = cfl * dx;
   int i;
   FILE *fp;

   // Grid and point values
   for(i=0; i<n; ++i)
   {
      x[i] = xmin + i*dx;
      up0[i] = uexact0(x[i]);
   }

   // Cell averages
   for(i=0; i<m; ++i)
      ua[i] = uexact1(x[i],x[i+1]);

   savesol(x, up0, ua);

   fp = fopen("sol0.dat","w");
   for(i=0; i<m; ++i)
      fprintf(fp,"%e %e\n", 0.5*(x[i]+x[i+1]), ua[i]);
   fclose(fp);

   for(i=0; i<m; ++i)
   {
      double um = midvalue(up0[i], up0[i+1], ua[i]);
      double xm = avgloc(up0[i], um, up0[i+1]);
      printf("%d %e\n", i, xm);
   }

   double t = 0.0;
   while(t < Tf)
   {
      if(t+dt > Tf) dt = Tf - t;

      // update to n+1/2
      for(i=0; i<n; ++i)
      {
         up1[i] = charint(i, x[i], up0, ua, 0.5*dt);
      }

      // update to n+1
      for(i=0; i<n; ++i)
      {
         up2[i] = charint(i, x[i], up0, ua, dt);
      }

      // Update cell average
      for(i=0; i<m; ++i)
      {
         double fl = (up0[i] + 4*up1[i] + up2[i]) / 6.0;
         double fr = (up0[i+1] + 4*up1[i+1] + up2[i+1]) / 6.0;
         ua[i] = ua[i] - (dt/dx) * ( fr - fl);
      }

      for(i=0; i<n; ++i)
         up0[i] = up2[i];

      t += dt;
      printf("t, dt = %e %e\n", t, dt);
   }

   fp = fopen("sol1.dat","w");
   for(i=0; i<m; ++i)
      fprintf(fp,"%e %e\n", 0.5*(x[i]+x[i+1]), ua[i]);
   fclose(fp);

   return 0;
}
