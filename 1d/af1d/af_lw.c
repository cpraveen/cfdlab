// Lax-Wendroff for face values
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define min(a,b)     ( (a) < (b) ? (a) : (b) )
#define max(a,b)     ( (a) > (b) ? (a) : (b) )
#define min3(a,b,c)  ( min(a,min(b,c)) )
#define sign(a)      ( (a) < 0 ? -1 : 1)

#include "sine.h"

#define NPOINTS 100

const int n = NPOINTS; // number of points
const int m = NPOINTS-1; // number of cells
const double xmin = -1.0, xmax = 1.0;
const double dx = (xmax - xmin)/(n-1);
const double Tf = 2.0;
const double cfl = 0.45;
const int tolimit = 0;
const double a = 1.0; // speed

// Simpson weights
const double w0 = 1.0/6.0, w1 = 4.0/6.0, w2 = 1.0/6.0;

double minmod(double a, double b, double c)
{
   int sa = sign(a);
   int sb = sign(b);
   int sc = sign(c);
   if(sa == sb && sb == sc) return sa*min3(fabs(a),fabs(b),fabs(c));
   return 0.0;
}

void savesol(double *x, double *up, double *um)
{
   int i, j;
   FILE *fp;

   fp = fopen("sol.dat", "w");
   for(i=0; i<m; ++i)
   {
      double xl = x[i];
      double xr = x[i+1];
      double xm = 0.5*(xl + xr);
      fprintf(fp,"%e %e\n", xl, up[i]);
      fprintf(fp,"%e %e\n", xm, um[i]);
      fprintf(fp,"%e %e\n", xr, up[i+1]);
   }
   fclose(fp);
}

// Unlimited scheme: compact stencil (-dx/2, +dx/2)
void flux1a(const double a, const double dt,
           const double *up0, const double *um, double *up1, double *fl)
{
   // update face values to n+1
   // left boundary
   {
      double ux = (um[0] - um[m-1])/dx;
      double uxx = (um[0] - 2.0*up0[0] + um[m-1])/(0.25*dx*dx);
      double ut = - a*ux;
      double utt = a*a*uxx;
      up1[0] = up0[0] + dt*ut + 0.5*dt*dt*utt;
      fl[0] = a*(up0[0] + 0.5*dt*ut + dt*dt*utt/6.0);
   }
   // interior faces
   for(int i=1; i<n-1; ++i)
   {
      double ux = (um[i] - um[i-1])/dx;
      double uxx = (um[i] - 2.0*up0[i] + um[i-1])/(0.25*dx*dx);
      double ut = - a*ux;
      double utt = a*a*uxx;
      up1[i] = up0[i] + dt*ut + 0.5*dt*dt*utt;
      fl[i] = a*(up0[i] + 0.5*dt*ut + dt*dt*utt/6.0);
   }
   // right boundary
   {
      up1[n-1] = up1[0];
      fl[n-1]  = fl[0];
   }
}

// Unlimited scheme: stencil (-dx,+dx)
void flux1b(const double a, const double dt,
           const double *up0, const double *um, double *up1, double *fl)
{
   // update face values to n+1
   // left boundary
   {
      double ux = (up0[1] - up0[n-2]) / (2*dx);
      double uxx = (up0[n-2] - 2.0 * up0[0] + up0[1]) / (dx * dx);
      double ut = -a * ux;
      double utt = a * a * uxx;
      up1[0] = up0[0] + dt * ut + 0.5 * dt * dt * utt;
      fl[0] = a * (up0[0] + 0.5 * dt * ut + dt * dt * utt / 6.0);
   }
   // interior faces
   for (int i = 1; i < n - 1; ++i)
   {
      double ux = (up0[i+1] - up0[i-1]) / (2*dx);
      double uxx = (up0[i-1] - 2.0 * up0[i] + up0[i+1]) / (dx * dx);
      double ut = -a * ux;
      double utt = a * a * uxx;
      up1[i] = up0[i] + dt * ut + 0.5 * dt * dt * utt;
      fl[i] = a * (up0[i] + 0.5 * dt * ut + dt * dt * utt / 6.0);
   }
   // right boundary
   {
      up1[n - 1] = up1[0];
      fl[n - 1] = fl[0];
   }
}

// Limited scheme
void flux2(const double a, const double dt,
           const double *up0, const double *um, double *up1, double *fl)
{
   const double beta = 1.0;

   // update face values to n+1
   // left boundary
   {
      double uxl = -(-1.5*up0[0] + 2.0*um[m-1] - 0.5*up0[n-2])/(0.5*dx);
      double uxc =  (um[0] - um[m-1])/dx;
      double uxr =  (-1.5*up0[0] + 2.0*um[0] - 0.5*up0[1])/(0.5*dx);
      double ux = minmod(uxc, beta*uxl, beta*uxr);

      double uxxl = (up0[n-2] - 2.0*um[m-1] + up0[n-1])/(0.25*dx*dx);
      double uxxc = (um[0] - 2.0*up0[0] + um[m-1])/(0.25*dx*dx);
      double uxxr = (up0[0] - 2.0*um[0] + up0[1])/(0.25*dx*dx);
      double uxx = minmod(uxxc, beta*uxxl, beta*uxxr);

      double ut = - a*ux;
      double utt = a*a*uxx;
      up1[0] = up0[0] + dt*ut + 0.5*dt*dt*utt;
      fl[0] = a*(up0[0] + 0.5*dt*ut + dt*dt*utt/6.0);
   }
   // interior faces
   for(int i=1; i<n-1; ++i)
   {
      double uxl = -(-1.5*up0[i] + 2.0*um[i-1] - 0.5*up0[i-1])/(0.5*dx);
      double uxc =  (um[i] - um[i-1])/dx;
      double uxr =  (-1.5*up0[i] + 2.0*um[i] - 0.5*up0[i+1])/(0.5*dx);
      double ux = minmod(uxc, beta*uxl, beta*uxr);

      double uxxl = (up0[i-1] - 2.0*um[i-1] + up0[i])/(0.25*dx*dx);
      double uxxc = (um[i] - 2.0*up0[i] + um[i-1])/(0.25*dx*dx);
      double uxxr = (up0[i] - 2.0*um[i] + up0[i+1])/(0.25*dx*dx);
      double uxx = minmod(uxxc, beta*uxxl, beta*uxxr);

      double ut = - a*ux;
      double utt = a*a*uxx;
      up1[i] = up0[i] + dt*ut + 0.5*dt*dt*utt;
      fl[i] = a*(up0[i] + 0.5*dt*ut + dt*dt*utt/6.0);
   }
   // right boundary
   {
      up1[n-1] = up1[0];
      fl[n-1]  = fl[0];
   }
}

int main()
{

   double x[n], up0[n], up1[n], um[n], ua[m], fl[n];

   double dt = cfl * dx / a;
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
   {
      double xm = 0.5*(x[i] + x[i+1]);
      um[i] = uexact0(xm);
      ua[i] = w0 * up0[i] + w1 * um[i] + w2 * up0[i+1];
   }

   savesol(x, up0, um);

   fp = fopen("sol0.dat","w");
   for(i=0; i<m; ++i)
      fprintf(fp,"%e %e\n", 0.5*(x[i]+x[i+1]), ua[i]);
   fclose(fp);

   double t = 0.0;
   while(t < Tf)
   {
      if(t+dt > Tf) dt = Tf - t;

      // update face values to n+1
      flux2(a, dt, up0, um, up1, fl);

      // Update cell average
      for(i=0; i<m; ++i)
      {
         ua[i] = ua[i] - (dt/dx) * ( fl[i+1] - fl[i]);
         um[i] = (ua[i] - w0*up1[i] - w2*up1[i+1])/w1;
      }

      for(i=0; i<n; ++i)
         up0[i] = up1[i];

      t += dt;
      printf("t, dt = %e %e\n", t, dt);
   }

   fp = fopen("sol1.dat","w");
   for(i=0; i<m; ++i)
      fprintf(fp,"%e %e\n", 0.5*(x[i]+x[i+1]), ua[i]);
   fclose(fp);

   savesol(x, up0, um);

   return 0;
}
