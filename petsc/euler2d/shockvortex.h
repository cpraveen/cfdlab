/* Shock vortex interaction
 * Run up to a time of 0.5 units
 *    ./fdweno -cfl 0.9 -da_grid_x 100 -da_grid_y 100 -Tf 0.5 -ts_monitor
 */

#define BC_LEFT    farfield
#define BC_RIGHT   farfield
#define BC_BOTTOM  neumann
#define BC_TOP     neumann

const double xmin = 0.0, xmax = 1.0;
const double ymin = 0.0, ymax = 1.0;
const double gas_gamma = 1.4;
const double gas_const = 1.0;
const int has_exact_sol = 0;
const double final_time = 0.5;

// Isentropic vortex
void exactsol(const double t, const double x1, const double y1, double *Prim)
{
   PetscPrintf(PETSC_COMM_WORLD,"exactsol not implemented\n");
   abort();
}

void initcond(const double x, const double y, double *Prim)
{
   const double x0   = 0.50;
   const double xc   = 0.25;
   const double yc   = 0.50;
   const double beta = 0.204;
   const double ep   = 0.3;
   const double rc   = 0.05;

   const double rhol = 1.0;
   const double ul   = sqrt(gas_gamma);
   const double vl   = 0.0;
   const double pl   = 1.0;

   const double pr   = 1.3;
   const double rhor = rhol*((gas_gamma-1.0)+(gas_gamma+1.0)*pr)/((gas_gamma+1.0)+(gas_gamma-1.0)*pr);
   const double ur   = sqrt(gas_gamma)+sqrt(2.0)*((1-pr)/sqrt(gas_gamma-1.0+pr*(gas_gamma+1.0)));
   const double vr   = 0.0;

   if(x <= x0)
   {
      const double r =  (pow(x-xc,2) + pow(y-yc,2))/pow(rc,2);
      const double du = ep*((y-yc)/rc)*exp(beta*(1.0-r));
      const double dv = -ep*((x-xc)/rc)*exp(beta*(1.0-r));
      const double dtheta = -((gas_gamma-1.0)/(4.0*beta*gas_gamma))*pow(ep,2)*exp(2.0*beta*(1.0-r));
      Prim[0] = pow(pl/rhol + dtheta, 1.0/(gas_gamma-1.0));
      Prim[1] = ul + du;
      Prim[2] = vl + dv;
      Prim[3] = pow(pl/rhol + dtheta, gas_gamma/(gas_gamma-1.0));
   }
   else
   {
      Prim[0] = rhor;
      Prim[1] = ur;
      Prim[2] = vr;
      Prim[3] = pr;
   }
}

void boundary_value(const double t, const double x, const double y, double *Con)
{
   double Prim[nvar];
   initcond(x, y, Prim);
   prim2con(Prim, Con);
}
