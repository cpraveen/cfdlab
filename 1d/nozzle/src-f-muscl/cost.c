#include<stdio.h>

int
main ()
{
   double c, d, cmin, cmax, dmin, dmax, dc, dd;
   double ain, aout;
   int nc, nd, i, j;
   FILE *fpt;

   ain  = 1.0512;
   aout = 1.75;
   cmin = 0.4;
   cmax = 1.2;
   dmin = 3.8;
   dmax = 4.2;
   nc = 51;
   nd = 51;
   dc = (cmax - cmin) / (nc-1);
   dd = (dmax - dmin) / (nd-1);

   system("rm -f cost");
   for (i = 0; i < nc; i++){
      for (j = 0; j < nd; j++) {
         c = cmin + dc * i;
         d = dmin + dd * j;
         printf("%d %d %lf %lf\n", i, j, c, d);
         fpt = fopen ("inp.dat", "w");
         fprintf (fpt, "%lf\n", ain);
         fprintf (fpt, "%lf\n", aout);
         fprintf (fpt, "%lf\n", c);
         fprintf (fpt, "%lf\n", d);
         fclose (fpt);
         system ("./nozzle > log");
         system ("tail -n 1 log >> cost");
      }
      system("echo '' >> cost");
   }

}
