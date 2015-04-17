#include<stdio.h>

__global__ void VecAdd(float* A, float* B, float *C)
{
   int idx = blockIdx.x*blockDim.x + threadIdx.x;
   C[idx] = A[idx] + B[idx];
}

int main()
{
   int N=16*16*16;
   float *A=0, *B=0, *C=0;
   float *Ad=0, *Bd=0, *Cd=0;

   A = (float*)malloc(N*sizeof(float));
   B = (float*)malloc(N*sizeof(float));
   C = (float*)malloc(N*sizeof(float));
   if(0==A||0==B||0==C){
      printf("Could not allocate host memory\n");
      return 1;
   }

   for(int i=0; i<N; i++)
   {
      A[i] = i;
      B[i] = i;
   }

   cudaMalloc( (void**)&Ad, N*sizeof(float) );
   cudaMalloc( (void**)&Bd, N*sizeof(float) );
   cudaMalloc( (void**)&Cd, N*sizeof(float) );

   if(0==Ad||0==Bd||0==Cd){
      printf("Could not allocate device memory\n");
      return 2;
   }

   cudaMemcpy(Ad, A, N*sizeof(float), cudaMemcpyHostToDevice);
   cudaMemcpy(Bd, B, N*sizeof(float), cudaMemcpyHostToDevice);

   dim3 grid, block;
   block.x = 1;
   grid.x  = N/block.x;

   VecAdd<<<grid, block>>>(Ad, Bd, Cd);

   cudaMemcpy(C, Cd, N*sizeof(float), cudaMemcpyDeviceToHost);

   for(int i=0; i<N; i++)
      printf("%d %e %e %e\n", i, A[i], B[i], C[i]);

   free(A); free(B); free(C);
   cudaFree(Ad); cudaFree(Bd); cudaFree(Cd);

   return 0;
}
