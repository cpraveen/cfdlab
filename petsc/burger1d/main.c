static char help[] = "Solves 1d Burger equation.\n\n";

#include <petscsys.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscvec.h>

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))

const double ark[] = {0.0, 3.0/4.0, 1.0/3.0};
const double xmin = 0.0;
const double xmax = 1.0;

double initcond(double x)
{
   return 1.0 + 0.5 * sin(2*M_PI*x);
}

//------------------------------------------------------------------------------
// Weno reconstruction of Jiang-Shu
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
double Flux(double u)
{
   return 0.5*pow(u,2);
}
//------------------------------------------------------------------------------
// Numerical flux
//------------------------------------------------------------------------------
double numflux(double ul, double ur)
{
   return max(Flux(max(0,ul)),Flux(min(0,ur)));
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
         fprintf(f, "%e %e\n", xmin+i*dx, uarray[i]);
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
   Vec      ug, ul;
   PetscInt i, ibeg, nloc, nx=100;
   const PetscInt sw = 3, ndof = 1; // stencil width
   PetscMPIInt rank, size;
   double cfl = 0.4;

   ierr = PetscInitialize(&argc, &argv, (char*)0, help); CHKERRQ(ierr);

   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
   MPI_Comm_size(PETSC_COMM_WORLD, &size);

   ierr = DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, -nx, ndof, sw, NULL, &da); CHKERRQ(ierr);

   ierr = DMCreateGlobalVector(da, &ug); CHKERRQ(ierr);

   ierr = DMDAGetCorners(da, &ibeg, 0, 0, &nloc, 0, 0); CHKERRQ(ierr);
   ierr = DMDAGetInfo(da,0,&nx,0,0,0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);
   PetscReal dx = (xmax - xmin) / (PetscReal)(nx);
   PetscPrintf(PETSC_COMM_WORLD,"nx = %d, dx = %e\n", nx, dx);
   PetscReal umax = 0;
   for(i=ibeg; i<ibeg+nloc; ++i)
   {
      PetscReal x = xmin + i*dx;
      PetscReal v = initcond(x);
      umax = max(umax,fabs(v));
      ierr = VecSetValues(ug,1,&i,&v,INSERT_VALUES); CHKERRQ(ierr);
   }
   ierr = VecAssemblyBegin(ug);  CHKERRQ(ierr);
   ierr = VecAssemblyEnd(ug);    CHKERRQ(ierr);

   savesol(nx, dx, ug);

   // Get local view
   ierr = DMGetLocalVector(da, &ul); CHKERRQ(ierr);

   PetscInt il, nl;
   ierr = DMDAGetGhostCorners(da,&il,0,0,&nl,0,0); CHKERRQ(ierr);

   double res[nloc], uold[nloc];
   double dt = cfl * dx / umax;
   double lam= dt/dx;

   double tfinal = 0.25, t = 0.0;

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

         if(rk==0)
            for(i=ibeg; i<ibeg+nloc; ++i) uold[i-ibeg] = u[i];

         for(i=0; i<nloc; ++i) 
            res[i] = 0.0;

         // Loop over faces and compute flux
         for(i=0; i<nloc+1; ++i)
         {
            // face between j-1, j
            int j   = il+sw+i;
            int jm1 = j-1;
            int jm2 = j-2;
            int jm3 = j-3;
            int jp1 = j+1;
            int jp2 = j+2;
            double ul = weno5(u[jm3],u[jm2],u[jm1],u[j],u[jp1]);
            double ur = weno5(u[jp2],u[jp1],u[j],u[jm1],u[jm2]);
            double flux = numflux(ul, ur);
            if(i==0)
            {
               res[i] -= flux;
            }
            else if(i==nloc)
            {
               res[i-1] += flux;
            }
            else
            {
               res[i]   -= flux;
               res[i-1] += flux;
            }
         }

         // Update solution
         for(i=ibeg; i<ibeg+nloc; ++i)
            unew[i] = ark[rk]*uold[i-ibeg] + (1-ark[rk])*(u[i] - lam * res[i-ibeg]);

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
