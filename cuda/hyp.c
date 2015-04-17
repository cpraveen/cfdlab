/* Linear convection equation with periodic BC
 * solved using MUSCL scheme
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define minmod(a,b)   ( (fabs(a) < fabs(b)) ? (a) : (b) )

//function prototypes
void initCond(int, float, float*);
void fluxFun(float*, float*);

float airk[3] = {0.0, 3.0/4.0, 1.0/3.0};

int main(){
   float *uold, *u, *fl, res;
   int np = 101, ns;
   float dx = 1.0/(np-1);
   float dt, cfl;
   int niter, maxiter;
   int nirk, rkmax=3;
   int i;
   FILE *fpt;

   ns = np + 2 + 2;

   uold = (float*)malloc(ns*sizeof(float));
   u    = (float*)malloc(ns*sizeof(float));
   fl   = (float*)malloc((np+1)*sizeof(float));

   cfl = 0.9;
   dt  = cfl*dx;
   maxiter = 1.0/dt + 1;

   //set initial conditions
   initCond(np, dx, u);

   fpt = fopen("init.dat", "w");
   for(i=0; i<np; i++) fprintf(fpt, "%e %e\n", dx*i, u[i+2]);
   fclose(fpt);


   //time-step loop
   for(niter=0; niter<maxiter; niter++){

      for(i=0; i<ns; i++) uold[i] = u[i];

      //RK stages
      for(nirk=0; nirk<rkmax; nirk++){

         //flux computation
         for(i=0; i<np+1; i++) fluxFun(&u[i+1], &fl[i]);

         //update conserved variable
         for(i=0; i<np; i++){
            res = fl[i+1] - fl[i];
            u[i+2] = airk[nirk]*uold[i+2] + 
                     (1.0-airk[nirk])*(u[i+2] - (dt/dx)*res);
         }

         //set periodicity
         u[0]    = u[np];
         u[1]    = u[np+1];
         u[np+2] = u[2];
         u[np+3] = u[3];

      }

   }

   fpt = fopen("final.dat", "w");
   for(i=0; i<np; i++) fprintf(fpt, "%e %e\n", dx*i, u[i+2]);
   fclose(fpt);

   free(uold);
   free(u);
   free(fl);

}

//set initial condition
void initCond(int np, float dx, float *u){
   int i;
   float x;

   for(i=0; i<np; i++){
      x = dx*i;
      u[i+2] = sin(2.0*M_PI*x);
   }

   //set periodicity
   u[0]    = u[np];
   u[1]    = u[np+1];
   u[np+2] = u[2];
   u[np+3] = u[3];

}

//flux function
void fluxFun(float *u, float *fl){
   float uj, ujp1, ujm1, ujp2;
   float ul, ur;

   uj   = *(u);
   ujp1 = *(u+1);
   ujm1 = *(u-1);
   ujp2 = *(u+2);

   ul = uj   + 0.5*minmod( (uj-ujm1), (ujp1-uj) );
   ur = ujp1 - 0.5*minmod( (ujp1-uj), (ujp2-ujp1) );

   *fl = ul;

}
