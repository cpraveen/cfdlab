/* Program to visualize a grid in amdba format using vigie. This reads the
   amdba file, creates a file in the format understood by vigie and starts
   vigie to visualize the grid.
 */
#include<stdio.h>
#include<stdlib.h>

/* --------------------------------------------------------------------- */
int main(int argc, char *argv[])
/* --------------------------------------------------------------------- */
{
   int i, j, dummy, FirstNodeNumber = 1, no_grid_points, ntri;
   int **tri;
   double *x, *y;
   FILE *fptr, *fvigie, *fdesc;

   if(argc < 2) {
      printf("Usage: %s <amdba file>\n", argv[0]);
      exit(0);
   }
   fptr = fopen(argv[1], "r");
   if(fptr == NULL) {
      printf("Could not open file %s\n", argv[1]);
      exit(0);
   }
   printf("Reading grid data in AMDBA format ...... ");
   fscanf(fptr, "%d%d", &no_grid_points, &ntri);
   x = (double *) calloc(no_grid_points, sizeof(double));
   y = (double *) calloc(no_grid_points, sizeof(double));
   if(x == NULL || y == NULL) {
      printf("Memory allocation was not successful\n");
      exit(0);
   }
   for(i = 0; i < no_grid_points; i++) {
      fscanf(fptr, "%d%lf%lf%d", &dummy, &x[i], &y[i], &dummy);
   }
   tri = (int **) calloc(ntri, sizeof(int *));
   for(i = 0; i < ntri; i++) {
      tri[i] = (int *) calloc(3, sizeof(int));
      fscanf(fptr, "%d%d%d%d%d", &dummy, &tri[i][0], &tri[i][1], &tri[i][2],
             &dummy);
      if(tri[i][0] == 0 || tri[i][1] == 0 || tri[i][2] == 0)
         FirstNodeNumber = 0;
   }
   fclose(fptr);

   if(FirstNodeNumber == 0)
      for(i = 0; i < ntri; i++) {
         tri[i][0] = tri[i][0] + 1;
         tri[i][1] = tri[i][1] + 1;
         tri[i][2] = tri[i][2] + 1;
         if(tri[i][0] < 1 || tri[i][1] < 1 || tri[i][2] < 1) {
            printf("Triangs numbering is negative at triangle %d.\n", i);
            exit(0);
         }
      }

   printf("Done.\n");
   printf("Number of grid points = %d\n", no_grid_points);
   printf("Number of triangles   = %d\n", ntri);

   /* Write grid in vigie format */
   fvigie = fopen("/tmp/vg.dat", "w");
   fprintf(fvigie, "points   %d\n", no_grid_points);
   for(i = 0; i < no_grid_points; i++)
      fprintf(fvigie, "%15.6f %15.6f\n", x[i], y[i]);
   fprintf(fvigie, "triangles    %d\n", ntri);
   for(i = 0; i < ntri; i++)
      fprintf(fvigie, "%8d %8d %8d\n", tri[i][0], tri[i][1], tri[i][2]);
   fprintf(fvigie, "end_block\n");
   fclose(fvigie);

   /* Create a description file for vigie */
   fdesc = fopen("/tmp/vg.desc", "w");
   fprintf(fdesc, "ascii2d\n");
   fprintf(fdesc, "/tmp/vg.dat\n");
   fclose(fdesc);

   /* Start vigie */
   system("vigie /tmp/vg.desc &");
}
