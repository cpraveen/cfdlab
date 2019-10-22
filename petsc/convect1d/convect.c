static char help[] = "Solves u_t + u_x = 0.\n\n";

#include <petscsys.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscvec.h>

const double ark[] = {0.0, 3.0/4.0, 1.0/3.0};
const double xmin = -1.0;
const double xmax = +1.0;

double initcond(double x)
{
   if(x >= -0.8 && x <= -0.6)
      return exp(-log(2)*pow(x+0.7,2)/0.0009);
   else if(x >= -0.4 && x <= -0.2)
      return 1.0;
   else if(x >= 0.0 && x <= 0.2)
      return 1.0 - fabs(10*(x-0.1));
   else if(x>= 0.4 && x <= 0.6)
      return sqrt(1 - 100*pow(x-0.5,2));
   else
      return 0.0;
}

//------------------------------------------------------------------------------
// Weno reconstruction
// Return left value for face between u0, up1
//------------------------------------------------------------------------------
double weno5(double um2, double um1, double u0, double up1, double up2)
{
   double eps = 1.0e-6;
   double gamma1=1.0/10.0, gamma2=3.0/5.0, gamma3=3.0/10.0;
   double beta1, beta2, beta3;
   double u1, u2, u3;
   double w1, w2, w3;

   beta1 = (13.0/12.0)*pow((um2 - 2.0*um1 + u0),2) +
           (1.0/4.0)*pow((um2 - 4.0*um1 + 3.0*u0),2);
   beta2 = (13.0/12.0)*pow((um1 - 2.0*u0 + up1),2) +
           (1.0/4.0)*pow((um1 - up1),2);
   beta3 = (13.0/12.0)*pow((u0 - 2.0*up1 + up2),2) +
           (1.0/4.0)*pow((3.0*u0 - 4.0*up1 + up2),2);

   w1 = gamma1 / pow(eps+beta1, 2);
   w2 = gamma2 / pow(eps+beta2, 2);
   w3 = gamma3 / pow(eps+beta3, 2);

   u1 = (1.0/3.0)*um2 - (7.0/6.0)*um1 + (11.0/6.0)*u0;
   u2 = -(1.0/6.0)*um1 + (5.0/6.0)*u0 + (1.0/3.0)*up1;
   u3 = (1.0/3.0)*u0 + (5.0/6.0)*up1 - (1.0/6.0)*up2;

   return (w1 * u1 + w2 * u2 + w3 * u3)/(w1 + w2 + w3);
}

//------------------------------------------------------------------------------
// Save solution to file
//------------------------------------------------------------------------------
PetscErrorCode savesol(int nx, double dx, Vec ug)
{
   int i, rank;
   static int count = 0;
   PetscErrorCode ierr;

   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   // Gather entire solution on rank=0 process. This is bad thing to do
   // in a real application.
   VecScatter ctx;
   Vec        uall;
   ierr = VecScatterCreateToZero(ug,&ctx,&uall); CHKERRQ(ierr);
   // scatter as many times as you need
   ierr = VecScatterBegin(ctx,ug,uall,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
   ierr = VecScatterEnd(ctx,ug,uall,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
   // destroy scatter context and local vector when no longer needed
   ierr = VecScatterDestroy(&ctx); CHKERRQ(ierr);
   if(rank==0)
   {
      PetscScalar *uarray;
      ierr = VecGetArray(uall, &uarray); CHKERRQ(ierr);
      FILE *f;
      f = fopen("sol.dat","w");
      for(i=0; i<nx; ++i)
         fprintf(f, "%e %e\n", xmin+i*dx+0.5*dx, uarray[i]);
      fclose(f);
      printf("Wrote solution into sol.dat\n");
      ierr = VecRestoreArray(uall, &uarray); CHKERRQ(ierr);
      if(count==0)
      {
         // Initial solution is copied to sol0.dat
         system("cp sol.dat sol0.dat");
         count = 1;
      }
   }
   ierr = VecDestroy(&uall); CHKERRQ(ierr);
   return(0);
}

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
   PetscErrorCode ierr;
   DM       da;
   Vec      ug; // global vector
   Vec      ul; // local vector
   PetscInt i, ibeg, nloc;
   PetscInt nx=200;         // no. of cells, can change via command line
   const PetscInt sw = 3;   // stencil width
   const PetscInt ndof = 1; // no. of dofs per cell
   PetscMPIInt rank, size;
   double cfl = 0.4;

   ierr = PetscInitialize(&argc, &argv, (char*)0, help); CHKERRQ(ierr);

   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
   MPI_Comm_size(PETSC_COMM_WORLD, &size);

   ierr = DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, nx, ndof, sw, NULL, &da); CHKERRQ(ierr);
   ierr = DMSetFromOptions(da); CHKERRQ(ierr);
   ierr = DMSetUp(da); CHKERRQ(ierr);
   ierr = DMCreateGlobalVector(da, &ug); CHKERRQ(ierr);

   ierr = DMDAGetCorners(da, &ibeg, 0, 0, &nloc, 0, 0); CHKERRQ(ierr);
   ierr = DMDAGetInfo(da,0,&nx,0,0,0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);
   PetscReal dx = (xmax - xmin) / (PetscReal)(nx);
   PetscPrintf(PETSC_COMM_WORLD,"nx = %d, dx = %e\n", nx, dx);
   for(i=ibeg; i<ibeg+nloc; ++i)
   {
      PetscReal x = xmin + i*dx + 0.5*dx;
      PetscReal v = initcond(x);
      ierr = VecSetValues(ug,1,&i,&v,INSERT_VALUES); CHKERRQ(ierr);
   }
   ierr = VecAssemblyBegin(ug);  CHKERRQ(ierr);
   ierr = VecAssemblyEnd(ug);    CHKERRQ(ierr);

   // save initial condition
   savesol(nx, dx, ug);

   // Get local view
   ierr = DMGetLocalVector(da, &ul); CHKERRQ(ierr);

   PetscInt il, nl;
   ierr = DMDAGetGhostCorners(da,&il,0,0,&nl,0,0); CHKERRQ(ierr);

   double res[nloc], uold[nloc];
   double dt = cfl * dx;
   double lam= dt/dx;

   double tfinal = 2.0, t = 0.0;

   while(t < tfinal)
   {
      if(t+dt > tfinal)
      {
         dt = tfinal - t;
         lam = dt/dx;
      }
      for(int rk=0; rk<3; ++rk)
      {
         ierr = DMGlobalToLocalBegin(da, ug, INSERT_VALUES, ul); CHKERRQ(ierr);
         ierr = DMGlobalToLocalEnd(da, ug, INSERT_VALUES, ul); CHKERRQ(ierr);

         PetscScalar *u;
         ierr = DMDAVecGetArrayRead(da, ul, &u); CHKERRQ(ierr);

         PetscScalar *unew;
         ierr = DMDAVecGetArray(da, ug, &unew); CHKERRQ(ierr);

         // First stage, store solution at time level n into uold
         if(rk==0)
            for(i=ibeg; i<ibeg+nloc; ++i) uold[i-ibeg] = u[i];

         for(i=0; i<nloc; ++i) 
            res[i] = 0.0;

         // Loop over faces and compute flux
         for(i=0; i<nloc+1; ++i) // local index
         {
            // face between j-1, j
            int j   = il+sw+i; // global index
            int jm1 = j-1;
            int jm2 = j-2;
            int jm3 = j-3;
            int jp1 = j+1;
            double uleft = weno5(u[jm3],u[jm2],u[jm1],u[j],u[jp1]);
            double flux = uleft;
            if(i==0) // first face
            {
               res[i] -= flux;
            }
            else if(i==nloc) // last face
            {
               res[i-1] += flux;
            }
            else // intermediate faces
            {
               res[i]   -= flux;
               res[i-1] += flux;
            }
         }

         // Update solution
         for(i=ibeg; i<ibeg+nloc; ++i)
            unew[i] = ark[rk]*uold[i-ibeg] + (1.0-ark[rk])*(u[i] - lam * res[i-ibeg]);

         ierr = DMDAVecRestoreArrayRead(da, ul, &u); CHKERRQ(ierr);
         ierr = DMDAVecRestoreArray(da, ug, &unew); CHKERRQ(ierr);
      }

      t += dt;
      PetscPrintf(PETSC_COMM_WORLD,"t = %f\n", t);
   }

   savesol(nx, dx, ug);

   // Destroy everything before finishing
   ierr = DMDestroy(&da); CHKERRQ(ierr);
   ierr = VecDestroy(&ug); CHKERRQ(ierr);

   ierr = PetscFinalize(); CHKERRQ(ierr);
}
