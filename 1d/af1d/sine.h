double uexact0(double x)
{
   return sin(2*M_PI*x);
}

double uexact1(double xl, double xr)
{
   return (cos(2*M_PI*xl) - cos(2*M_PI*xr))/(2*M_PI*(xr-xl));
}
