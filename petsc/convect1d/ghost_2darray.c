#include<stdlib.h>
#include<stdio.h>

// No of unknowns at each grid point/cell
#define nvar  3

int main()
{
   int nx = 10; // no of real cells
   int ng = 3;  // no of ghost cells on either side
   int n  = ng + nx + ng; // total no of cells
   double (*v)[nvar] = calloc(n, sizeof(*v)) + ng;

   // Real cells
   for(int i=0; i<nx; ++i)
      for(int j=0; j<nvar; ++j)
         v[i][j] = i + j;

   // Fill with periodicity
   // fill ghost values at left boundary
   for(int i=-ng; i<0; ++i)
      for(int j=0; j<nvar; ++j)
         v[i][j] = v[i+nx][j];
   // fill ghost values at right boundary
   for(int i=nx; i<nx+ng; ++i)
      for(int j=0; j<nvar; ++j)
         v[i][j] = v[i-nx][j];

   // Print all cells
   for(int i=-ng; i<nx+ng; ++i)
      printf("%d  %f %f %f\n", i, v[i][0], v[i][1], v[i][2]);
   return 0;
}
