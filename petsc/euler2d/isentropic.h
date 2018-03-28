#define BC_LEFT   periodic
#define BC_RIGHT  periodic
#define BC_BOTTOM periodic
#define BC_TOP    periodic

const double xmin = -5.0, xmax = 5.0;
const double ymin = -5.0, ymax = 5.0;
const double gas_gamma = 1.4;
const double gas_const = 1.0;

// Isentropic vortex
void exactsol(const double t, const double x1, const double y1, double *Prim)
{
   const double M = 0.5;
   const double alpha = 0.0; // between 0 and 90
   const double beta = 5.0;
   const double vx0 = M * cos(alpha*M_PI/180);
   const double vy0 = M * sin(alpha*M_PI/180);

   double x = x1 - vx0 * t;
   double y = y1 - vy0 * t;

   while(x < -5.0) x += 10.0;
   while(y < -5.0) y += 10.0;

   const double r2 = x*x + y*y;
   Prim[0] =  pow(1.0 - (gas_gamma-1.0)*(beta*beta)/(8.0*gas_gamma*M_PI*M_PI)*exp(1-r2), (1.0/(gas_gamma-1.0)));
   Prim[1] =  vx0 - beta/(2.0*M_PI)*y*exp(0.5*(1.0-r2));
   Prim[2] =  vy0 + beta/(2.0*M_PI)*x*exp(0.5*(1.0-r2));
   Prim[3] =  pow(Prim[0],gas_gamma);
}

void initcond(const double x, const double y, double *Prim)
{
   exactsol(0.0, x, y, Prim);
}

void boundary_value(const double t, const double x, const double y, double *Con)
{
   abort();
}
