#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double smooth(double x)
{
   return sin(2*M_PI*x);
}

double smooth_x(double x)
{
   return 2*M_PI*cos(2*M_PI*x);
}

double sine4(double x)
{
   return pow(sin(2*M_PI*x),4);
}

double sine4_x(double x)
{
   return 4*pow(sin(2*M_PI*x),3)*2*M_PI*cos(2*M_PI*x);
}

double hat(double x)
{
   if(x < 0.25 || x > 0.75)
      return 0.0;
   else
      return 1.0;
}

double hat_x(double x)
{
   return 0.0;
}

int main()
{
   const double a = 1.0;
   const double Tf = 10.0;
   const double cfl = 0.95;
   const int n = 101;
   const int m = n-1;
   const double xmin = 0.0, xmax = 1.0;
   const double dx = (xmax-xmin)/(n-1);
   double x[n], u[n][2], v[m][2];
   int i;
   FILE *fpt;

   // set initial condition
   for(i=0; i<n; ++i)
   {
      x[i] = xmin + i*dx;
      u[i][0] = sine4(x[i]);
      u[i][1] = 0.25*dx*sine4_x(x[i]);
   }

   const double dt = cfl * dx / fabs(a);
   const double nu = a * dt / dx;
   double t = 0.0;
   int it = 0;
   while(t < Tf)
   {
      // update to mid time level
      for(i=0; i<m; ++i)
      {
         v[i][0] = 0.5*(
                    (1+nu)*u[i][0] + (1-nu)*u[i+1][0] + 
                    (1-nu*nu)*(u[i][1] - u[i+1][1]));
         v[i][1] = -0.5*( u[i][0] - u[i+1][0] + 
                         (1-nu)*u[i][1] + (1+nu)*u[i+1][1]);
      }

      // update to next time level
      i = 0;
      u[i][0] = 0.5*(
                 (1+nu)*v[m-1][0] + (1-nu)*v[i][0] + 
                 (1-nu*nu)*(v[m-1][1] - v[i][1]));
      u[i][1] = -0.5*( v[m-1][0] - v[i][0] + 
                      (1-nu)*v[m-1][1] + (1+nu)*v[i][1]);
      for(i=1; i<n-1; ++i)
      {
         u[i][0] = 0.5*(
                    (1+nu)*v[i-1][0] + (1-nu)*v[i][0] + 
                    (1-nu*nu)*(v[i-1][1] - v[i][1]));
         u[i][1] = -0.5*( v[i-1][0] - v[i][0] + 
                         (1-nu)*v[i-1][1] + (1+nu)*v[i][1]);
      }
      i = n-1;
      u[i][0] = 0.5*(
                 (1+nu)*v[i-1][0] + (1-nu)*v[0][0] + 
                 (1-nu*nu)*(v[i-1][1] - v[0][1]));
      u[i][1] = -0.5*( v[i-1][0] - v[0][0] + 
                      (1-nu)*v[i-1][1] + (1+nu)*v[0][1]);

      t += dt; ++it;
      printf("%d  %e\n",it,t);
   }

   fpt = fopen("sol.dat","w");
   for(i=0; i<n; ++i)
   {
      double x1 = x[i] - a*t;
      double ue = sine4(x1);
      fprintf(fpt,"%e %e %e %e\n",x[i],u[i][0],ue,u[i][0]-ue);
   }
   fclose(fpt);
}
