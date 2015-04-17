#define DEFINE_GLOBALS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "decl.h"
#include "param.h"

#define NFACES 101
#define NCELLS (NFACES-1)

void residu(int, int, double *, double **, double **);
void costfunc(double *, double, double *);

int main()
{
   int iter, i, j, nc, nf;
   double x[NFACES], a[NFACES], **q, **qold, **res,
      dt, residue, residue1, maxres, ar, cost, dum, ptarg[NCELLS];
   FILE *fpt;

   nf = NFACES;
   nc = nf - 1;

   q = (double **) calloc(nc, sizeof(double *));
   qold = (double **) calloc(nc, sizeof(double *));
   res = (double **) calloc(nc, sizeof(double *));
   for(i = 0; i < nc; i++) {
      q[i] = (double *) calloc(3, sizeof(double));
      qold[i] = (double *) calloc(3, sizeof(double));
      res[i] = (double *) calloc(3, sizeof(double));
   }

   init();
   printf("Flux for flow solver = %d\n", flux1);
   read_shape(nf, x, a);
// init_shape(nf, x, a)

//Initialize flow
   for(i = 0; i < nc; i++) {
      q[i][0] = rinf;
      q[i][1] = rinf * uinf;
      q[i][2] = pinf / gam1 + 0.5 * rinf * uinf * uinf;
   }

   maxres = 1.0e-10;
   residue = residue1 = 1.0;
   iter = 0;
   fpt = fopen("res.dat", "w");

//Time step loop
   while(iter < maxiter && residue > maxres) {
      time_step(nc, q, &dt);
      save_old(nc, q, qold);

      residu(nc, nf, a, q, res);

      for(i = 0; i < nc; i++) {
         ar = 0.5 * (a[i] + a[i + 1]);
         for(j = 0; j < 3; j++)
            q[i][j] = qold[i][j] - dt * res[i][j] / dx / ar;
      }

      sol_residue(nc, res, &residue);
      iter = iter + 1;
      if(iter == 1)
         residue1 = residue;
      residue = residue / residue1;
      fprintf(fpt, "%d %e %e\n", iter, log10(residue), dt);
   }
   fclose(fpt);

   printf("Number of flow iterations = %d\n", iter);

//Read target pressure from a file
   fpt = fopen("flow-target.dat", "r");
   if(fpt==NULL){
      printf("File flow-target.dat not found !!!\n");
      exit(0);
   }
   for(i = 0; i < nc; i++)
      fscanf(fpt, "%lf%lf%lf%lf%lf", &dum, &dum, &dum, &ptarg[i], &dum);
   fclose(fpt);

//Compute cost function - L2 norm of pressure difference
   cost = 0.0;
   for(i = 0; i < nc; i++)
      costfunc(q[i], ptarg[i], &cost);
   printf("Cost function = %e\n", cost);
   fpt = fopen("cost.dat", "w");
   fprintf(fpt, "%e", cost);;
   fclose(fpt);
//    print*, dc, dd, cost

//Save solution into a file
   result(nf, nc, x, q);

   for(i = 0; i < nc; i++) {
      free(q[i]);
      free(qold[i]);
      free(res[i]);
   }
   free(q);
   free(qold);
   free(res);

   return 1;
}
