#define DEFINE_GLOBALS
#include "adolc.h"
#include <stdio.h>
#include <stdlib.h>
#include "decl.h"
#include "param.h"

#define NFACES 101
#define NCELLS (NFACES-1)

//using namespace std;

void residu(int, int, double *, adouble **, adouble **);
void residu(int, int, adouble *, double **, adouble **);
void costfunc(adouble *, double, adouble *);
void print_tapestats(int);

int main()
{
   int iter, ii, i, j, nc, nf;
   double x[NFACES], a[NFACES], **q, res_dummy,
      dt, residue, residue1, maxres, ar, dum, ptarg[NCELLS];
   FILE *fpt;
   short int tag;

   nf = NFACES;
   nc = nf - 1;

   q = (double **) calloc(nc, sizeof(double *));
   for(i = 0; i < nc; i++) {
      q[i] = (double *) calloc(3, sizeof(double));
   }

   double *costb = new double[1];
   double *cost_q = new double[3 * nc];

   adouble *bq[nc], *bres[nc];
   for(i = 0; i < nc; i++) {
      bq[i] = new adouble[3];
      bres[i] = new adouble[3];
   }
   adouble *ba = new adouble[nf];

   double *resb = new double[3 * nc];
   double *resb_old = new double[3 * nc];
   double *qb = new double[3 * nc];

   init();
   printf("Flux for flow solver = %d\n", flux1);
   read_shape(nf, x, a);
// init_shape(nf, x, a)

   read_flow(nc, q);

//Read target pressure from a file
   fpt = fopen("flow-target.dat", "r");
   if(fpt == NULL) {
      printf("File flow-target.dat not found !!!\n");
      exit(0);
   }
   for(i = 0; i < nc; i++)
      fscanf(fpt, "%lf%lf%lf%lf%lf", &dum, &dum, &dum, &ptarg[i], &dum);
   fclose(fpt);

//Gradient of cost function wrt flow
   tag = 0;
   trace_on(tag, 1);
   for(i = 0; i < nc; i++)
      for(j = 0; j < 3; j++)
         bq[i][j] <<= q[i][j];
   adouble cost = 0.0;
   for(i = 0; i < nc; i++)
      costfunc(bq[i], ptarg[i], &cost);
   cost >>= dum;
   trace_off();
   print_tapestats(tag);
   printf("Cost function = %e\n", cost.value());
   costb[0] = 1.0;
   reverse(tag, 1, 3 * nc, 0, costb, cost_q);

//Generate tape for adjoint equation
   tag = 0;
   trace_on(tag, 1);
   for(i = 0; i < nc; i++)
      for(j = 0; j < 3; j++)
         bq[i][j] <<= q[i][j];
   residu(nc, nf, a, bq, bres);
   for(i = 0; i < nc; i++)
      for(j = 0; j < 3; j++)
         bres[i][j] >>= res_dummy;
   trace_off();
   print_tapestats(tag);

   maxres = 1.0e-10;
   residue = residue1 = 1.0;
   iter = 0;
   fpt = fopen("res.dat", "w");

//Time step loop
   for(i = 0; i < 3 * nc; i++)
      resb[i] = 0.0;
   time_step(nc, q, &dt);
   while(iter < maxiter && residue > maxres) {
      for(i = 0; i < 3 * nc; i++)
         resb_old[i] = resb[i];

      //reverse(tag, 3 * nc, 3 * nc, 0, resb, qb);
      fos_reverse(tag, 3 * nc, 3 * nc, resb, qb);

      for(i = 0; i < nc; i++) {
         ar = 0.5 * (a[i] + a[i + 1]);
         for(j = 0; j < 3; j++) {
            ii = 3 * i + j;
            qb[ii] = qb[ii] + cost_q[ii];
            resb[ii] = resb_old[ii] - dt * qb[ii] / dx / ar;
         }
      }

      residue = 0.0;
      for(i = 0; i < 3 * nc; i++)
         residue += qb[i] * qb[i];
      residue = sqrt(residue) / nc;
      iter = iter + 1;
      fprintf(fpt, "%d %e\n", iter, log10(residue));
   }
   fclose(fpt);

   printf("Number of adjoint iterations = %d\n", iter);
   fpt = fopen("flowb.dat", "w");
   for(i = 0; i < nc; i++)
      fprintf(fpt, "%18.10e%18.10e%18.10e%18.10e\n", x[i] + 0.5 * dx,
              resb[3 * i], resb[3 * i + 1], resb[3 * i + 2]);
   fclose(fpt);

//Now we compute gradient trans(adj) * dR/da
//Generate tape for adjoint equation
   //for(i = 0; i < nf; i++)
   //ba[i] = a[i];
   tag = 0;
   trace_on(tag, 1);
   for(i = 1; i < nf - 1; i++)
      ba[i] <<= a[i];
   residu(nc, nf, ba, q, bres);
   for(i = 0; i < nc; i++)
      for(j = 0; j < 3; j++)
         bres[i][j] >>= res_dummy;
   trace_off();
   printf("Tape for dR/da\n");
   print_tapestats(tag);

   fos_reverse(tag, 3 * nc, nf - 2, resb, qb);

   //Save gradient
   fpt = fopen("shapeb.dat", "w");
   fprintf(fpt, "%18.10e%18.10e\n", x[0], 0.0);
   for(i = 1; i < nf - 1; i++)
      fprintf(fpt, "%18.10e%18.10e\n", x[i], qb[i - 1]);
   fprintf(fpt, "%18.10e%18.10e\n", x[nf - 1], 0.0);
   fclose(fpt);

   for(i = 0; i < nc; i++) {
      free(q[i]);
   }
   free(q);

   return 1;
}

/* Print tape size statistics */
void print_tapestats(int tag)
{
   int count[100];
   tapestats(tag, count);
   printf("Number of independents         = %6d\n", count[0]);
   printf("Number of dependents           = %6d\n", count[1]);
   printf("Max number of live active vars = %6d\n", count[2]);
   printf("Size of value stack            = %6d\n", count[3]);
   printf("Buffer size                    = %6d\n", count[4]);
}
