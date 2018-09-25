// Double Mach reflection problem

#define BC_LEFT    farfield
#define BC_RIGHT   neumann
#define BC_BOTTOM  wall
#define BC_TOP     farfield 

const double xmin = 0.0, xmax = 4.0;
const double ymin = 0.0, ymax = 1.0;
const double gas_gamma = 1.4;
const double gas_const = 1.0;
const int has_exact_sol = 0;
const double final_time = 0.2;

void exactsol(const double t, const double x1, const double y1, double *Prim)
{
   PetscPrintf(PETSC_COMM_WORLD,"exactsol not implemented\n");
   abort();
}

// This must be same as boundary_value for t=0.
void initcond(const double x, const double y, double *Prim)
{
   const double x_0 = 1.0/6.0;
   const double sq3 = 1.0/sqrt(3.0);

   if(x < x_0 + y*sq3)
   {
      Prim[0] = 8.0;
      Prim[1] = 4.125/sq3;// 8.25cos(pi/6)
      Prim[2] = -4.125;//    8.25sin(pi/6)
      Prim[3] = 116.5;
   }
   else
   {
      Prim[0] = 1.4;
      Prim[1] = 0.0;
      Prim[2] = 0.0;
      Prim[3] = 1.0;
   }
}

// NOTE: This must return conserved variables.
// This must be used only on the top boundary.
void boundary_value(const double t, const double x, const double y, double *Con)
{
   // ensure that bc for t=0 also holds...
   double shock_pos;
   const double x_0 = 1.0/6.0;

   shock_pos = x_0 + (y + 20.0*t)/sqrt(3.0);
   
   double Prim[nvar];

   if(x < shock_pos)
   {
      Prim[0] = 8.0;
      Prim[1] = 4.125*sqrt(3.0);
      Prim[2] = -4.125;
      Prim[3] = 116.5;
   }
   else
   {
      Prim[0] = 1.4;
      Prim[1] = 0.0;
      Prim[2] = 0.0;
      Prim[3] = 1.0;
   }

   prim2con(Prim, Con);
}
