#define BC_LEFT   neumann
#define BC_RIGHT  neumann
#define BC_BOTTOM neumann
#define BC_TOP    neumann

const double xmin = 0.0, xmax = 1.0;
const double ymin = 0.0, ymax = 1.0;
const double gas_gamma = 1.4;
const double gas_const = 1.0;
const int has_exact_sol = 0;
const double final_time = 0.8;

// 2-D Riemann
void exactsol(const double t, const double x1, const double y1, double *Prim)
{
   PetscPrintf(PETSC_COMM_WORLD,"exactsol not implemented\n");
   abort();
}

void initcond(const double x, const double y, double *Prim)
{
   if((x >= 0.0 && x <= 1.0) && (y >= 0.0 && y <= 1.0))
   {
      if(x >= 0.8 && y >= 0.8)
      {
         Prim[0] = 1.5;
         Prim[1] = 0.0;
         Prim[2] = 0.0;
         Prim[3] = 1.5;
      }
      else if(x < 0.8 && y >= 0.8)
      {
         Prim[0] = 0.5323;
         Prim[1] = 1.206;
         Prim[2] = 0.0;
         Prim[3] = 0.3;
      }
      else if(x < 0.8 && y < 0.8)
      {
         Prim[0] = 0.138;
         Prim[1] = 1.206;
         Prim[2] = 1.206;
         Prim[3] = 0.029;
      }
      else
      {
         Prim[0] = 0.5323;
         Prim[1] = 0.0;
         Prim[2] = 1.206;
         Prim[3] = 0.3;
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
