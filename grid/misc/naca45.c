/* NACA 4 and 5 series airfoils */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#define npmax 1000
#define CHORD 1.00893041136514
#define RADIUS 15.0

#define SOLID 3
#define OUTER 5

int npts, opts;
double m, p, t, c;
double xr, yr, area;

int main(int argc, char *argv[])
{
   int np, series, code;
   double xp[npmax], yc[npmax], dyc[npmax], yt[npmax], x[npmax], y[npmax],
      xo[npmax], yo[npmax];
   void xcoord();
   void naca4series();
   void naca5series();
   void thickness(int, double[], double[]);
   void total();
   void writepoints(double[], double[]);
   void outerpoints(double[], double[]);
   void writepoly(double[], double[], double[], double[]);
   void writepts(double[], double[], double[], double[]);

   if(argc < 4) {
      printf("Usage  : %s <NACA airfoil> <Int1> <Int2>\n", argv[0]);
      printf("         Int1 = Number of points on airfoil\n");
      printf("         Int2 = Number of points on outer boundary\n");
      printf("Example: %s 2415 100 30\n", argv[0]);
      exit(0);
   }

   series = strlen(argv[1]);
   code = atoi(argv[1]);
   np = atoi(argv[2]);
   opts = atoi(argv[3]);

   xcoord(np, xp);

   switch (series) {
      case 4:
         m = atoi(&argv[1][0]) / 1000;
         p = atoi(&argv[1][1]) / 100;
         t = atoi(&argv[1][2]);
         naca4series(np, xp, yc, dyc);
         break;
      case 5:
         c = atoi(&argv[1][0]) / 10000;
         p = atoi(&argv[1][1]) / 100;
         t = atoi(&argv[1][3]);
         naca5series(np, xp, yc, dyc);
         break;
      default:
         printf("Only NACA 4 and 5 series airfoils supported.\n");
         exit(0);
   }


   thickness(np, xp, yt);
   total(np, xp, yc, dyc, yt, x, y);
   writepoints(x, y);

   outerpoints(xo, yo);
   xr = 0.5;
   yr = yc[np / 2];
   writepts(x, y, xo, yo);
   writepoly(x, y, xo, yo);

}

void xcoord(int np, double xp[])
{
   int i;
   double theta;

   for(i = 0; i < np; i++) {
      theta = M_PI * i / (np - 1);
      xp[i] = 0.5 * (1.0 + cos(theta));
   }

}

/* NACA 4 series airfoil */
void naca4series(int np, double x[], double yc[], double dyc[])
{
   int i;
   m = m / 100.0;
   p = p / 10.0;
   t = t / 100.0;
   printf("NACA 4 series airfoil\n");
   printf("Maximum camber         = %lf\n", m);
   printf("Position of max camber = %lf\n", p);
   printf("Maximum thickness      = %lf\n", t);

   if(m == 0.0) {
      for(i = 0; i < np; i++) {
         yc[i] = 0.0;
         dyc[i] = 0.0;
      }
      return;
   }

   for(i = 0; i < np; i++)
      if(x[i] <= p) {
         yc[i] = m * (2.0 * p * x[i] - x[i] * x[i]) / p / p;
         dyc[i] = 2.0 * m * (p - x[i]) / p / p;
      }
      else {
         yc[i] =
            m * ((1.0 - 2.0 * p) + 2.0 * p * x[i] -
                 x[i] * x[i]) / (1.0 - p) / (1.0 - p);
         dyc[i] = 2.0 * m * (p - x[i]) / (1.0 - p) / (1.0 - p);
      }


}

/* NACA 5 series airfoil */
void naca5series(int np, double x[], double yc[], double dyc[])
{
   int i;
   double m, k1;

   c = 1.5 * c / 10.0;
   p = p / 2.0 / 100.0;
   t = t / 100.0;

   if(p == 0.05) {
      m = 0.0580;
      k1 = 361.4;
   }
   else if(p == 0.10) {
      m = 0.1260;
      k1 = 51.64;
   }
   else if(p == 0.15) {
      m = 0.2025;
      k1 = 15.957;
   }
   else if(p == 0.20) {
      m = 0.2900;
      k1 = 6.643;
   }
   else if(p == 0.25) {
      m = 0.3910;
      k1 = 3.23;
   }
   else {
      printf("Error: Naca 5 series\n");
      printf("Dont know how to find m and k1 for p = %f\n", &p);
      printf("Known airfoils are 210xx, 220xx, 230xx, 240xx and 250xx\n");
      exit(0);
   }

   printf("NACA 5 series airfoil\n");
   printf("c = %lf\n", c);
   printf("m = %lf\n", m);
   printf("p = %lf\n", p);
   printf("t = %lf\n", t);
   printf("k1= %lf\n", k1);

   for(i = 0; i < np; i++)
      if(x[i] <= p) {
         yc[i] =
            (k1 / 6.0) * (x[i] * x[i] * x[i] - 3.0 * m * x[i] * x[i] +
                          m * m * (3.0 - m) * x[i]);
         dyc[i] =
            (k1 / 6.0) * (3.0 * x[i] * x[i] - 6.0 * m * x[i] +
                          m * m * (3.0 - m));
      }
      else {
         yc[i] = k1 * m * m * m * (1.0 - x[i]) / 6.0;
         dyc[i] = -k1 * m * m * m / 6.0;
      }
}

/* Thickness distribution */
void thickness(int np, double x[], double yt[])
{
   int i;
   double c1, c2, c3, c4, c5, xx;

   c1 = 0.2969;
   c2 = 0.1260;
   c3 = 0.3516;
   c4 = 0.2843;
   c5 = 0.1015;
   for(i = 0; i < np; i++) {
      xx = CHORD * x[i];
      yt[i] =
         5.0 * t * (c1 * sqrt(xx) - c2 * xx - c3 * xx * xx +
                    c4 * xx * xx * xx - c5 * xx * xx * xx * xx);
   }
}

/* Add thickness and camber to get final profile */
void
total(int np, double xp[], double yc[], double dyc[], double yt[],
      double x[], double y[])
{
   int i;
   double theta;

   npts = 0;

   /* Lower surface */
   for(i = 0; i < np; i++) {
      theta = atan(dyc[i]);
      x[npts] = xp[i] + yt[i] * sin(theta);
      y[npts] = yc[i] - yt[i] * cos(theta);
      ++npts;
   }

   /* Upper surface */
   for(i = np - 2; i > 0; i--) {
      theta = atan(dyc[i]);
      x[npts] = xp[i] - yt[i] * sin(theta);
      y[npts] = yc[i] + yt[i] * cos(theta);
      ++npts;
   }


   /* Trailing edge is not at y=0; shift it */
   for(i = 1; i < npts; i++)
      y[i] = y[i] - y[0];
   y[0] = 0.0;

   for(i = 2; i < npts - 1; i++) {
      x[i - 1] = x[i];
      y[i - 1] = y[i];
   }
   x[npts - 2] = x[0];
   y[npts - 2] = y[0];
   npts = npts - 1;

   printf("Total number of points = %d\n", npts);

}

/* Write profile to a file */
void writepoints(double x[], double y[])
{
   int i;
   FILE *fpt;

   fpt = fopen("naca.dat", "w");
   for(i = 0; i < npts; i++)
      fprintf(fpt, "%20.10e %20.10e\n", x[i], y[i]);
   fclose(fpt);
}

/* Outer boundary points */
void outerpoints(double xo[], double yo[])
{
   int i;
   double dtheta, theta, dx, dy, lent;

   dtheta = 2.0 * M_PI / opts;

   for(i = 0; i < opts; i++) {
      theta = dtheta * i;
      xo[i] = RADIUS * cos(theta);
      yo[i] = RADIUS * sin(theta);
   }

   dx = xo[0] - xo[1];
   dy = yo[0] - yo[1];
   lent = dx * dx + dy * dy;
   area = 0.25 * sqrt(3.0) * lent;
}

/* Poly file for Shewchuk's triangle */
void writepoly(double x[], double y[], double xo[], double yo[])
{
   int i, nseg;
   FILE *fpoly;

   printf("Writing points in .poly format -> naca.poly\n");
   printf("Triangle usage:\n");
   printf("\ttriangle -pq30a%f naca.poly\n", area);

   fpoly = fopen("naca.poly", "w");
   fprintf(fpoly, "%d %d %d %d\n", npts - 1 + opts, 2, 0, 1);
   for(i = 0; i < npts - 1; i++)
      fprintf(fpoly, "%5d %16.8e %16.8e %5d\n", i + 1, x[i], y[i], SOLID);
   for(i = 0; i < opts; i++)
      fprintf(fpoly, "%5d %16.8e %16.8e %5d\n", i + npts, xo[i], yo[i],
              OUTER);
   nseg = npts - 1 + opts;
   fprintf(fpoly, "%5d %5d\n", nseg, 1);
   for(i = 0; i < npts - 2; i++)
      fprintf(fpoly, "%5d %5d %5d %5d\n", i + 1, i + 1, i + 2, 3);
   fprintf(fpoly, "%5d %5d %5d\n", npts - 1, npts - 1, 1, 3);
   for(i = 0; i < opts - 1; i++)
      fprintf(fpoly, "%5d %5d %5d %5d\n", i + npts, i + npts, i + npts + 1,
              5);
   fprintf(fpoly, "%5d %5d %5d %5d\n", npts + opts - 1, npts + opts - 1,
           npts, 5);
   fprintf(fpoly, "%5d\n", 1);
   fprintf(fpoly, "%5d %16.8e %16.8e\n", 1, xr, yr);
   fclose(fpoly);
}

/* pts file for delaundo */
void writepts(double x[], double y[], double xo[], double yo[])
{
   int i;
   FILE *fpts;

   printf("Writing points in .dpl format  -> naca.pts\n");

   fpts = fopen("naca.pts", "w");
   fprintf(fpts, "NEWBND\n");
   fprintf(fpts, "NAMEBN\n");
   fprintf(fpts, "1\n");
   fprintf(fpts, "NFRSBN\n");
   fprintf(fpts, "1\n");
   fprintf(fpts, "NLSTBN\n");
   fprintf(fpts, "1\n");
   fprintf(fpts, "ITYPBN\n");
   fprintf(fpts, "1\n");
   fprintf(fpts, "NRBNDE\n");
   fprintf(fpts, "%d\n", npts);
   fprintf(fpts, "BNDEXY\n");
   for(i = 0; i < npts; i++)
      fprintf(fpts, "%18.8e %18.8e\n", x[i], y[i]);
   fprintf(fpts, "NEWBND\n");
   fprintf(fpts, "NAMEBN\n");
   fprintf(fpts, "2\n");
   fprintf(fpts, "NFRSBN\n");
   fprintf(fpts, "2\n");
   fprintf(fpts, "NLSTBN\n");
   fprintf(fpts, "2\n");
   fprintf(fpts, "ITYPBN\n");
   fprintf(fpts, "2\n");
   fprintf(fpts, "NRBNDE\n");
   fprintf(fpts, "%d\n", opts + 1);
   fprintf(fpts, "BNDEXY\n");
   for(i = 0; i < opts; i++)
      fprintf(fpts, "%18.8e %18.8e\n", xo[i], yo[i]);
   fprintf(fpts, "%18.8e %18.8e\n", xo[0], yo[0]);
   fprintf(fpts, "ENDDAT");
   fclose(fpts);

}
