// Authors: Aman Saxena, Praveen C
#define BC_LEFT   periodic
#define BC_RIGHT  periodic
#define BC_BOTTOM periodic
#define BC_TOP    periodic

const double xmin = 0.0, xmax = 1.0;
const double ymin = 0.0, ymax = 1.0;
const double gas_gamma = 1.4;
const double gas_const = 1.0;
const int has_exact_sol = 0;
const double final_time = 0.8;

//2-D Riemann
void exactsol(const double t, const double x1, const double y1, double *Prim)
{
   PetscPrintf(PETSC_COMM_WORLD,"exactsol not implemented\n");
   abort();
}

void initcond(const double x, const double y, double *Prim)
{
   double w, sigma, z_1, z_2;

   w     = 0.1;
   sigma = 0.05/sqrt(2.0);
   z_1   = -0.5*pow((y-0.25)/sigma,2);
   z_2   = -0.5*pow((y-0.75)/sigma,2);

   if((x >= 0.0 && x <= 1.0 )&& (y >= 0.0 && y <= 1.0))
   {
      Prim[3] = 2.5;
      Prim[2] = w*sin(4.0*M_PI*x)*(exp(z_1) + exp(z_2));

      if(y <= 0.75 && y > 0.25)
      {
         Prim[0] = 2.0;
         Prim[1] = 0.5;
      }
      else
      {
         Prim[0] = 1.0;
         Prim[1] = -0.5;
      }
   }
   else
   {
      PetscPrintf(PETSC_COMM_WORLD,"initcond undefined for x,y\n");
      abort();
   }
}

void boundary_value(const double t, const double x, const double y, double *Con)
{
   PetscPrintf(PETSC_COMM_WORLD,"boundary_value not implemented\n");
   abort();
}
