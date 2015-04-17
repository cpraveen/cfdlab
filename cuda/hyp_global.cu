/* Linear convection equation with periodic BC
 * solved using MUSCL scheme
 * CUDA implementation of hyp.c using only global memory
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define minmod(a,b)   ( (fabs(a) < fabs(b)) ? (a) : (b) )


//function prototypes
void initCond(int, float, float*);
__global__ void fluxFun(float*, float*);
__global__ void update(int, float, float, float*, float*, float*);
__global__ void periodic(float*, int);

int main(){
   float *u;
   float *uold_d, *u_d, *fl_d;
   int np = 101, ns;
   float dx = 1.0/(np-1);
   float dt, cfl;
   int niter, maxiter;
   int nirk, rkmax=3;
   int i;
   FILE *fpt;
   dim3 grid, block;

   ns   = np + 2 + 2;

   u    = (float*)malloc(ns*sizeof(float));

   cudaMalloc((void**)&uold_d, (ns)*sizeof(float));
   cudaMalloc((void**)&u_d,    (ns)*sizeof(float));
   cudaMalloc((void**)&fl_d,   (np+1)*sizeof(float));

   cfl = 0.9;
   dt  = cfl*dx;
   maxiter = 1.0/dt + 1;

   //set initial conditions
   initCond(np, dx, u);

   fpt = fopen("init.dat", "w");
   for(i=0; i<np; i++) fprintf(fpt, "%e %e\n", dx*i, u[i+2]);
   fclose(fpt);

   cudaMemcpy(uold_d, u, (ns)*sizeof(float), 
              cudaMemcpyHostToDevice);
   cudaMemcpy(u_d, u, (ns)*sizeof(float), 
              cudaMemcpyHostToDevice);

   //time-step loop
   for(niter=0; niter<maxiter; niter++){

      //RK stages
      for(nirk=0; nirk<rkmax; nirk++){

         //flux computation
         block.x = 3;
         grid.x  = (np+1)/block.x;
         fluxFun<<<grid,block>>>(u_d, fl_d);

         //update conserved variable
         block.x = 1;
         grid.x  = (np)/block.x;
         update<<<grid,block>>>(nirk, dt, dx, uold_d, fl_d, u_d);

         //set periodicity
         block.x = 1;
         grid.x  = (ns)/block.x;
         periodic<<<grid,block>>>(u_d, np);
      }
      cudaMemcpy(uold_d, u_d, (ns)*sizeof(float), 
                 cudaMemcpyDeviceToDevice);

   }
   cudaMemcpy(u, u_d, (ns)*sizeof(float), 
              cudaMemcpyDeviceToHost);

   fpt = fopen("final.dat", "w");
   for(i=0; i<np; i++) fprintf(fpt, "%e %e\n", dx*i, u[i+2]);
   fclose(fpt);

   free(u);

   cudaFree(uold_d);
   cudaFree(u_d);
   cudaFree(fl_d);

}

//set initial condition
void initCond(int np, float dx, float *u){
   int i;
   float x;

   for(i=0; i<np; i++){
      x = dx*i;
      u[i+2] = sin(2.0*M_PI*x);
   }
   u[0]    = u[np];
   u[1]    = u[np+1];
   u[np+2] = u[2];
   u[np+3] = u[3];

}

//flux function
__global__ void fluxFun(float *u, float *fl){
   float uj, ujp1, ujm1, ujp2;
   float ul, ur;

   int idx = blockIdx.x*blockDim.x + threadIdx.x;

   ujm1 = *(u+idx);
   uj   = *(u+idx+1);
   ujp1 = *(u+idx+2);
   ujp2 = *(u+idx+3);

   ul = uj   + 0.5*minmod( (uj-ujm1), (ujp1-uj) );
   ur = ujp1 - 0.5*minmod( (ujp1-uj), (ujp2-ujp1) );

   fl[idx] = ul;

}

//perform one stage of RK
__global__ void update(int nirk, float dt, float dx, float *uold, float *fl, 
                       float *u){
   int idx = blockIdx.x*blockDim.x + threadIdx.x;
   float res;
   float airk[3] = {0.0, 3.0/4.0, 1.0/3.0};

   res = fl[idx+1] - fl[idx];
   u[idx+2] = airk[nirk]*uold[idx+2] + 
          (1.0-airk[nirk])*(u[idx+2] - (dt/dx)*res);
}

//set periodic BC
__global__ void periodic(float *u, int np){
   int idx = blockIdx.x*blockDim.x + threadIdx.x;

   if(idx==0)
      u[idx] = u[np];
   else if(idx==1)
      u[idx] = u[np+1];
   else if(idx==np+2)
      u[idx] = u[2];
   else if(idx==np+3)
      u[idx] = u[3];
   else
      u[idx] = u[idx];

}
