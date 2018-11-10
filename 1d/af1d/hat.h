double uexact0(double x)
{
   if(x > 0.25 && x < 0.75)
      return 1.0;
   else
      return 0.0;
}

double uexact1(double xl, double xr)
{
   if(xl < 0.25 && xr > 0.25)
   {
      return (xr - 0.25)/(xr - xl);
   }
   else if(xl < 0.75 && xr > 0.75)
   {
      return (0.75 - xl)/(xr - xl);
   }
   else
      return uexact0(0.5*(xl+xr));
}
