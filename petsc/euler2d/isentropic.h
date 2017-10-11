#define PERIODIC  2

const double xmin = -5.0, xmax = 5.0;
const double ymin = -5.0, ymax = 5.0;
const double gas_gamma = 1.4;
const double gas_const = 1.0;

// Isentropic vortex
void initcond(const double x, const double y, double *Prim)
{
   const double M = 0.5;
   const double alpha = 0.0;
   const double beta = 5.0;
   const double r2 = x*x + y*y;
   Prim[0] =  pow(1.0 - (gas_gamma-1.0)*(beta*beta)/(8.0*gas_gamma*M_PI*M_PI)*exp(1-r2), (1.0/(gas_gamma-1.0)));
   Prim[1] =  M*cos(alpha*M_PI/180.0) - beta/(2.0*M_PI)*y*exp(0.5*(1.0-r2));
   Prim[2] =  M*sin(alpha*M_PI/180.0) + beta/(2.0*M_PI)*x*exp(0.5*(1.0-r2));
   Prim[3] =  pow(Prim[0],gas_gamma);
}
