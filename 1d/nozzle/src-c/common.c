#include <math.h>
#include <stdio.h>
#include "param.h"
#include "decl.h"

//Initialize some constants and grid
void init()
{
   FILE *fpt;

   gam = 1.4;
   gam1 = gam - 1.0;
   gasconst = 1.0;

   minf = 1.5;
   rinf = 1.0;
   pinf = 1.0;
   Tinf = pinf / (rinf * gasconst);
   uinf = minf * sqrt(gam * gasconst * Tinf);

   prat = 2.5;
   pout = prat * pinf;

   printf("Inflow mach number = %f\n", minf);
   printf("Pressure ratio     = %f\n", prat);

   printf("Reading paramters from param.dat\n");
   fpt = fopen("param.dat", "r");
   fscanf(fpt, "%lf", &cfl);
   fscanf(fpt, "%d", &maxiter);
   fscanf(fpt, "%d%d", &flux1, &flux2);
   fclose(fpt);
}

//Copy array q into qold
void save_old(int nc, double **q, double **qold)
{
   int i, j;

   for(i = 0; i < nc; i++)
      for(j = 0; j < 3; j++)
         qold[i][j] = q[i][j];
}

//Convert conservative to primitive variables
void con_to_prim(double *r, double *u, double *p, double *q)
{
   *r = q[0];
   *u = q[1] / q[0];
   *p = gam1 * (q[2] - 0.5 * q[1] * q[1] / q[0]);
}

//Global time step
void time_step(int nc, double **q, double *dt)
{
   int i;
   double dt1, dtlocal, a, r, u, p;

   dt1 = 1.0;
   for(i = 0; i < nc; i++) {
      con_to_prim(&r, &u, &p, q[i]);
      a = sqrt(gam * p / r);
      dtlocal = dx / (fabs(u) + a);
      dt1 = (dtlocal < dt1) ? dtlocal : dt1;
   }
   dt1 = cfl * dt1;
   *dt = dt1;
}

//L2 norm of solution residual
void sol_residue(int nc, double **res, double *residue)
{
   int i;
   double res1;

   res1 = 0.0;
   for(i = 0; i < nc; i++)
      res1 +=
         res[i][0] * res[i][0] + res[i][1] * res[i][1] +
         res[i][2] * res[i][2];
   res1 = sqrt(res1) / nc;

   *residue = res1;
}

//Save flow solution and nozzle shape
void result(int nf, int nc, double *x, double **q)
{
   int i;
   double r, u, p, son, mach;
   FILE *fpt;

   fpt = fopen("flow.dat", "w");
   for(i = 0; i < nc; i++) {
      con_to_prim(&r, &u, &p, q[i]);
      son = sqrt(gam * p / r);
      mach = u / son;
      fprintf(fpt, "%18.10e%18.10e%18.10e%18.10e%18.10e\n", x[i] + 0.5 * dx,
              r, u, p, mach);
   }
   fclose(fpt);
}

//Read flow solution
void read_flow(int nc, double **q)
{
   int i;
   double r, u, p, dum, mach;
   FILE *fpt;

   printf("Reading flow solution from flow.dat\n");
   fpt = fopen("flow.dat", "r");
   for(i = 0; i < nc; i++) {
      fscanf(fpt, "%lf%lf%lf%lf%lf", &dum, &r, &u, &p, &mach);
      q[i][0] = r;
      q[i][1] = r * u;
      q[i][2] = p / gam1 + 0.5 * r * u * u;
   }
   fclose(fpt);
}

//Read nozzle shape
void read_shape(int nf, double *x, double *a)
{
   int i;
   FILE *fpt;

   printf("Reading nozzle shape from shape.dat\n");
   fpt = fopen("shape.dat", "r");
   for(i = 0; i < nf; i++)
      fscanf(fpt, "%lf%lf", &x[i], &a[i]);
   fclose(fpt);

   dx = x[1] - x[0];
   ain = a[0];
   aout = a[nf - 1];

   printf("Inlet area         = %e\n", ain);
   printf("Outlet area        = %e\n", aout);
   printf("Area ratio         = %e\n", aout / ain);
}

//Save adjoint flow solution and shape gradient
void result_adj(int nf, int nc, double *x, double *a, double **q)
{
   int i;
   FILE *fpt;

   fpt = fopen("flowb.dat", "w");
   for(i = 0; i < nc; i++)
      fprintf(fpt, "%18.10e%18.10e%18.10e%18.10e\n", x[i] + 0.5 * dx, q[i][0],
              q[i][1], q[i][2]);
   fclose(fpt);

   fpt = fopen("shapeb.dat", "w");
   for(i = 0; i < nf; i++)
      fprintf(fpt, "%18.10e%18.10e\n", x[i], a[i]);
   fclose(fpt);
}
