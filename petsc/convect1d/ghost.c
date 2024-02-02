#include<stdlib.h>
#include<stdio.h>

int main()
{
   int nx = 10; // no of real cells
   int ng = 3;  // no of ghost cells on each side
   int n  = ng + nx + ng; // total no of cells
   double *u = (double*)malloc(n * sizeof(double*)) + ng;

   // Real cells
   for(int i=0; i<nx; ++i)
      u[i] = i;

   // Fill with periodicity
   // fill ghost values at left boundary
   for(int i=-ng; i<0; ++i)
      u[i] = u[i+nx];
   // fill ghost values at right boundary
   for(int i=nx; i<nx+ng; ++i)
      u[i] = u[i-nx];

   // Print all cells
   for(int i=-ng; i<nx+ng; ++i)
      printf("%d  %e\n", i, u[i]);
   return 0;
}
