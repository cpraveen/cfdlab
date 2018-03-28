#define BC_LEFT    farfield
#define BC_RIGHT   supersonic
#define BC_BOTTOM  wall
#define BC_TOP     farfield

const double xmin = -5.0, xmax = 5.0;
const double ymin = -5.0, ymax = 5.0;
const double gas_gamma = 1.4;
const double gas_const = 1.0;

// Isentropic vortex
void exactsol(const double t, const double x1, const double y1, double *Prim)
{
   abort();
}

void initcond(const double x, const double y, double *Prim)
{
   Prim[0] = 1.0;
   Prim[1] = 2.9;
   Prim[2] = 0.0;
   Prim[3] = 1.0/gas_gamma;
}

void boundary_value(const double t, const double x, const double y, double *Con)
{
   Con[0] =  1.69997;
   Con[1] =  2.61934;
   Con[2] = -0.50633;
   Con[3] =  1.52819;
}
