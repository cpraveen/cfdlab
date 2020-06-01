/* Displays points from a *.pts file of delaundo using gnuplot
 * Use gnuplot 4.0 or later to allow zooming by mouse 
 * Written by Praveen. C, http://pc.freeshell.org */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define NPMAX  2000

#define NEWBND 0
#define NAMEBN 1
#define NFRSBN 2
#define NLSTBN 3
#define ITYPBN 4
#define NRBNDE 5
#define NRINDE 6
#define MINRES 7
#define BNDEXY 8
#define INDEXY 9
#define ENDDAT 10

/* ----------------------------------------------------------------------- */
int main(int argc, char *argv[])
/* ----------------------------------------------------------------------- */
{
   char *ptsfile, str[32];
   int i, stype, namebn, itypbn, nrbnde, nfrsbn, nlstbn, minres;
   int status;
   double x[NPMAX], y[NPMAX];
   FILE *fpt, *fpts, *fgnu;

   if(argc < 2) {
      printf("Usage: %s <pts file>\n", argv[0]);
      exit(0);
   }
   ptsfile = argv[1];
   fpt = fopen(ptsfile, "r");
   if(fpt == NULL) {
      printf("Could not open %s\n", ptsfile);
      exit(0);
   }
   fpts = fopen("/tmp/pts.dat", "w");

   while(1) {
      status = fscanf(fpt, "%s", str);
      if(status == EOF) {
         strcpy(str, "ENDDAT");
      }
      stype = type(str);
      switch (stype) {
         case (NEWBND):
            break;
         case (NAMEBN):
            fscanf(fpt, "%d", &namebn);
            break;
         case (NFRSBN):
            fscanf(fpt, "%d", &nfrsbn);
            break;
         case (NLSTBN):
            fscanf(fpt, "%d", &nlstbn);
            break;
         case (ITYPBN):
            fscanf(fpt, "%d", &itypbn);
            break;
         case (NRBNDE):
            fscanf(fpt, "%d", &nrbnde);
            break;
         case (NRINDE):
            fscanf(fpt, "%d", &nrbnde);
            break;
         case (MINRES):
            fscanf(fpt, "%d", &minres);
            break;
         case (BNDEXY):
         case (INDEXY):
            for(i = 0; i < nrbnde; i++) {
               fscanf(fpt, "%lf%lf", &x[i], &y[i]);
               fprintf(fpts, "%lf %lf\n", x[i], y[i]);
            }
            fprintf(fpts, "\n");
            printf("Boundary = %d, Type = %d, Points = %d\n",
                   namebn, itypbn, nrbnde);
            break;
         case (ENDDAT):
            fclose(fpt);
            fclose(fpts);
            fgnu = fopen("/tmp/pts.gnu", "w");
            fprintf(fgnu, "set size ratio -1\n");
            fprintf(fgnu, "set nokey\n");
            fprintf(fgnu, "plot '/tmp/pts.dat' w lp pt 6\n");
            fprintf(fgnu, "pause -1");
            fclose(fgnu);
            system("gnuplot /tmp/pts.gnu");
            exit(0);
      }
   }


}

int type(char str[])
{
   if(!strcmp(str, "NEWBND"))
      return NEWBND;
   if(!strcmp(str, "NAMEBN"))
      return NAMEBN;
   if(!strcmp(str, "NFRSBN"))
      return NFRSBN;
   if(!strcmp(str, "NLSTBN"))
      return NLSTBN;
   if(!strcmp(str, "ITYPBN"))
      return ITYPBN;
   if(!strcmp(str, "NRBNDE"))
      return NRBNDE;
   if(!strcmp(str, "NRINDE"))
      return NRINDE;
   if(!strcmp(str, "MINRES"))
      return MINRES;
   if(!strcmp(str, "BNDEXY"))
      return BNDEXY;
   if(!strcmp(str, "INDEXY"))
      return INDEXY;
   if(!strcmp(str, "ENDDAT"))
      return ENDDAT;
   printf("Fatal error: Unknown string %s\n", str);
   exit(0);
}
