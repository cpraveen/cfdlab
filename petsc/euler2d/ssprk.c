static char help[] = "Solves 2d Euler equations.\n\n";

#include <petscsys.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscvec.h>

#define min(a,b)  ( (a < b) ? a : b )
#define nvar  4

const PetscReal ark[] = {0.0, 3.0/4.0, 1.0/3.0};
const PetscReal xmin = -5.0, xmax = 5.0;
const PetscReal ymin = -5.0, ymax = 5.0;
const PetscReal gas_gamma = 1.4;
const PetscReal gas_const = 1.0;
PetscReal dx, dy;

// Isentropic vortex
void initcond(const PetscReal x, const PetscReal y, PetscReal *Prim)
{
   const PetscReal M = 0.5;
   const PetscReal alpha = 0.0;
   const PetscReal beta = 5.0;
   const PetscReal r2 = x*x + y*y;
   Prim[0] =  pow(1.0 - (gas_gamma-1.0)*(beta*beta)/(8.0*gas_gamma*M_PI*M_PI)*exp(1-r2), (1.0/(gas_gamma-1.0)));
   Prim[1] =  M*cos(alpha*M_PI/180.0) - beta/(2.0*M_PI)*y*exp(0.5*(1.0-r2));
   Prim[2] =  M*sin(alpha*M_PI/180.0) + beta/(2.0*M_PI)*x*exp(0.5*(1.0-r2));
   Prim[3] =  pow(Prim[0],gas_gamma);
}

//------------------------------------------------------------------------------
// Weno reconstruction
// Return left value for face between u0, up1
//------------------------------------------------------------------------------
PetscReal weno5(const PetscReal um2, const PetscReal um1, const PetscReal u0, 
             const PetscReal up1, const PetscReal up2)
{
   PetscReal eps = 1.0e-6;
   PetscReal gamma1=1.0/10.0, gamma2=3.0/5.0, gamma3=3.0/10.0;
   PetscReal beta1, beta2, beta3;
   PetscReal u1, u2, u3;
   PetscReal w1, w2, w3;

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

// Conserved to primitive variables
void con2prim(const PetscReal *Con, PetscReal *Prim)
{
  Prim[0] = Con[0];
  Prim[1] = Con[1]/Con[0];
  Prim[2] = Con[2]/Con[0];
  Prim[3] = (Con[3] - 0.5*Prim[0]*(pow(Prim[1],2) + pow(Prim[2],2)))*(gas_gamma-1.0);
}

// Primitive to conserved variables
void prim2con(const PetscReal *Prim, PetscReal *Con)
{
  Con[0] = Prim[0];
  Con[1] = Prim[0]*Prim[1];
  Con[2] = Prim[0]*Prim[2];
  Con[3] = 0.5*Prim[0]*(pow(Prim[1],2) + pow(Prim[2],2)) + Prim[3]/(gas_gamma-1.0);
}

// Compute maximum eigenvalue in direction (nx,ny)
PetscReal maxeigval(const PetscReal *Con, const PetscReal nx, const PetscReal ny)
{
   PetscReal Prim[nvar];
   con2prim(Con, Prim);
   const PetscReal u = fabs(Prim[1]*nx + Prim[2]*ny);
   const PetscReal a = sqrt(gas_gamma*Prim[3]/Prim[0]);
   return u + a;
}

// Compute local timestep
PetscReal dt_local(const PetscReal *Con)
{
   PetscReal Prim[nvar];
   con2prim(Con, Prim);
   const PetscReal a = sqrt(gas_gamma*Prim[3]/Prim[0]);
   const PetscReal sx = fabs(Prim[1]) + a;
   const PetscReal sy = fabs(Prim[2]) + a;
   return 1.0/(sx/dx + sy/dy);
}

// Simple average flux
void avg_flux(const PetscReal *Ul, const PetscReal *Ur, 
              const PetscReal nx, const PetscReal ny, PetscReal *flux)
{
   PetscReal Pl[nvar], Pr[nvar];
   con2prim(Ul, Pl);
   con2prim(Ur, Pr);

   PetscReal fluxl[nvar], fluxr[nvar];

  fluxl[0] = Ul[1]*nx + Ul[2]*ny;
  fluxl[1] = Pl[3]*nx + Pl[1]*fluxl[0];
  fluxl[2] = Pl[3]*ny + Pl[2]*fluxl[0];
  fluxl[3] = (Ul[3]+Pl[3])*(Pl[1]*nx + Pl[2]*ny);

  fluxr[0] = Ur[1]*nx + Ur[2]*ny;
  fluxr[1] = Pr[3]*nx + Pr[1]*fluxr[0];
  fluxr[2] = Pr[3]*ny + Pr[2]*fluxr[0];
  fluxr[3] = (Ur[3]+Pr[3])*(Pr[1]*nx + Pr[2]*ny);

  for(int i=0; i<nvar; ++i) flux[i] = 0.5*(fluxl[i] + fluxr[i]);
}

// Rusanov flux
void numflux(const PetscReal *Ul, const PetscReal *Ur, 
             const PetscReal nx, const PetscReal ny, PetscReal *flux)
{
   PetscReal Ua[nvar];
   for(int i=0; i<nvar; ++i) Ua[i] = 0.5*(Ul[i] + Ur[i]);
   PetscReal lam = maxeigval( Ua, nx, ny );
   avg_flux(Ul,Ur,nx,ny,flux);
   for(int i=0; i<nvar; ++i) flux[i] -= 0.5*lam*(Ur[i] - Ul[i]);
}

//------------------------------------------------------------------------------
PetscErrorCode savesol(int *c, PetscReal t, DM da, Vec ug)
{
   PetscErrorCode ierr;
   char           filename[32] = "sol";
   PetscMPIInt    rank;
   PetscInt       i, j, nx, ny, ibeg, jbeg, nlocx, nlocy;
   FILE           *fp;
   Vec            ul;
   PetscScalar    ***u;

   ierr = DMGetLocalVector(da, &ul); CHKERRQ(ierr);
   ierr = DMGlobalToLocalBegin(da, ug, INSERT_VALUES, ul); CHKERRQ(ierr);
   ierr = DMGlobalToLocalEnd(da, ug, INSERT_VALUES, ul); CHKERRQ(ierr);
   ierr = DMDAVecGetArrayDOFRead(da, ul, &u); CHKERRQ(ierr);
   ierr = DMDAGetInfo(da,0,&nx,&ny,0,0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);
   ierr = DMDAGetCorners(da, &ibeg, &jbeg, 0, &nlocx, &nlocy, 0); CHKERRQ(ierr);

   PetscInt iend = PetscMin(ibeg+nlocx+1, nx);
   PetscInt jend = PetscMin(jbeg+nlocy+1, ny);

   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
   sprintf(filename, "sol-%03d-%03d.plt", *c, rank);
   fp = fopen(filename,"w");
   fprintf(fp, "TITLE = \"u_t + u_x + u_y = 0\"\n");
   fprintf(fp, "VARIABLES = x, y, rho, u, v, p\n");
   fprintf(fp, "ZONE STRANDID=1, SOLUTIONTIME=%e, I=%d, J=%d, DATAPACKING=POINT\n", 
           t, iend-ibeg, jend-jbeg);
   for(j=jbeg; j<jend; ++j)
      for(i=ibeg; i<iend; ++i)
   {
      PetscReal x = xmin + i*dx + 0.5*dx;
      PetscReal y = ymin + j*dy + 0.5*dy;
      PetscReal prim[nvar];
      con2prim(u[j][i], prim);
      fprintf(fp, "%e %e %e %e %e %e\n", x, y, prim[0], prim[1], prim[2], prim[3]);
   }
   fclose(fp);

   ierr = DMDAVecRestoreArrayDOFRead(da, ul, &u); CHKERRQ(ierr);
   ierr = DMRestoreLocalVector(da, &ul); CHKERRQ(ierr);

   ++(*c);
   return(0);
}
//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
   // some parameters that can overwritten from command line
   PetscReal Tf  = 10.0;
   PetscReal cfl = 0.8;
   PetscInt  si  = 100;
   PetscInt  nx  = 50, ny=50; // use -da_grid_x, -da_grid_y to override these
   
   PetscErrorCode ierr;
   DM       da;
   Vec      ug, ul;
   PetscInt i, j, ibeg, jbeg, nlocx, nlocy, d;
   const PetscInt sw = 3; // stencil width
   PetscMPIInt rank, size;
   PetscScalar ***u;
   PetscScalar ***unew;
   int c = 0; // counter for saving solution files

   ierr = PetscInitialize(&argc, &argv, (char*)0, help); CHKERRQ(ierr);

   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
   MPI_Comm_size(PETSC_COMM_WORLD, &size);
   
   // Get some command line options
   ierr = PetscOptionsGetReal(NULL,NULL,"-Tf",&Tf,NULL); CHKERRQ(ierr);
   ierr = PetscOptionsGetReal(NULL,NULL,"-cfl",&cfl,NULL); CHKERRQ(ierr);
   ierr = PetscOptionsGetInt(NULL,NULL,"-si",&si,NULL); CHKERRQ(ierr);

   ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC,
                       DMDA_STENCIL_BOX, nx, ny, PETSC_DECIDE, PETSC_DECIDE, nvar,
                       sw, NULL, NULL, &da); CHKERRQ(ierr);
   ierr = DMSetFromOptions(da); CHKERRQ(ierr);
   ierr = DMSetUp(da); CHKERRQ(ierr);
   ierr = DMDAGetInfo(da,0,&nx,&ny,0,0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);
   dx = (xmax - xmin) / (PetscReal)(nx);
   dy = (ymax - ymin) / (PetscReal)(ny);
   PetscPrintf(PETSC_COMM_WORLD,"nx = %d, dx = %e\n", nx, dx);
   PetscPrintf(PETSC_COMM_WORLD,"ny = %d, dy = %e\n", ny, dy);

   ierr = DMCreateGlobalVector(da, &ug); CHKERRQ(ierr);
   ierr = PetscObjectSetName((PetscObject) ug, "Solution"); CHKERRQ(ierr);

   ierr = DMDAGetCorners(da, &ibeg, &jbeg, 0, &nlocx, &nlocy, 0); CHKERRQ(ierr);
   ierr = DMDAVecGetArrayDOF(da, ug, &u); CHKERRQ(ierr);
   for(j=jbeg; j<jbeg+nlocy; ++j)
      for(i=ibeg; i<ibeg+nlocx; ++i)
      {
         PetscReal x = xmin + i*dx + 0.5*dx;
         PetscReal y = ymin + j*dy + 0.5*dy;
         PetscReal prim[nvar];
         initcond(x, y, prim);
         prim2con(prim, u[j][i]);
      }
   ierr = DMDAVecRestoreArrayDOF(da, ug, &u); CHKERRQ(ierr);
   ierr = savesol(&c, 0.0, da, ug); CHKERRQ(ierr);

   // Get local view
   ierr = DMGetLocalVector(da, &ul); CHKERRQ(ierr);

   PetscInt il, jl, nl, ml;
   ierr = DMDAGetGhostCorners(da,&il,&jl,0,&nl,&ml,0); CHKERRQ(ierr);

   // Allocate res[nlocy][nlocx][nvar] and uold[nlocy][nlocx][nvar]
   PetscReal (*res) [nlocx][nvar] = calloc(nlocy, sizeof(*res) );
   PetscReal (*uold)[nlocx][nvar] = calloc(nlocy, sizeof(*uold));
   PetscReal dt, lam;
   PetscReal t = 0.0;
   int it = 0;

   while(t < Tf)
   {
      for(int rk=0; rk<3; ++rk)
      {
         // start global to local but continue with other tasks
         ierr = DMGlobalToLocalBegin(da, ug, INSERT_VALUES, ul); CHKERRQ(ierr);
         ierr = DMDAVecGetArrayDOF(da, ug, &unew); CHKERRQ(ierr);

         if(rk==0)
         {
            // compute time step
            PetscReal dtlocal = 1.0e20;

            for(j=jbeg; j<jbeg+nlocy; ++j)
               for(i=ibeg; i<ibeg+nlocx; ++i)
               {
                  for(d=0; d<nvar; ++d)
                     uold[j-jbeg][i-ibeg][d] = unew[j][i][d];
                  dtlocal = min(dtlocal, dt_local(unew[j][i]));
               }

            MPI_Allreduce(&dtlocal, &dt, 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);
            dt *= cfl;

            if(t+dt > Tf) dt = Tf - t;
            lam = dt/(dx*dy);
         }

         for(j=0; j<nlocy; ++j)
            for(i=0; i<nlocx; ++i)
               for(d=0; d<nvar; ++d)
                  res[j][i][d] = 0.0;

         // finish global to local
         ierr = DMGlobalToLocalEnd(da, ug, INSERT_VALUES, ul); CHKERRQ(ierr);
         ierr = DMDAVecGetArrayDOFRead(da, ul, &u); CHKERRQ(ierr);

         // x fluxes
         for(i=0; i<nlocx+1; ++i)
            for(j=0; j<nlocy; ++j)
            {
               // face between k-1, k
               PetscInt k = il+sw+i;
               PetscInt l = jl+sw+j;
               PetscReal UL[nvar], UR[nvar], flux[nvar];
               for(d=0; d<nvar; ++d)
               {
                  UL[d] = weno5(u[l][k-3][d],u[l][k-2][d],u[l][k-1][d],u[l][k][d],u[l][k+1][d]);
                  UR[d] = weno5(u[l][k+2][d],u[l][k+1][d],u[l][k][d],u[l][k-1][d],u[l][k-2][d]);
               }
               numflux(UL, UR, 1.0, 0.0, flux);
               if(i==0)
               {
                  for(d=0; d<nvar; ++d)
                     res[j][i][d] -= dy * flux[d];
               }
               else if(i==nlocx)
               {
                  for(d=0; d<nvar; ++d)
                     res[j][i-1][d] += dy * flux[d];
               }
               else
               {
                  for(d=0; d<nvar; ++d)
                  {
                     res[j][i][d]   -= dy * flux[d];
                     res[j][i-1][d] += dy * flux[d];
                  }
               }
            }

         // y fluxes
         for(j=0; j<nlocy+1; ++j)
            for(i=0; i<nlocx; ++i)
            {
               // face between l-1, l
               PetscInt k = il+sw+i;
               PetscInt l = jl+sw+j;
               PetscReal UL[nvar], UR[nvar], flux[nvar];
               for(d=0; d<nvar; ++d)
               {
                  UL[d] = weno5(u[l-3][k][d],u[l-2][k][d],u[l-1][k][d],u[l][k][d],u[l+1][k][d]);
                  UR[d] = weno5(u[l+2][k][d],u[l+1][k][d],u[l][k][d],u[l-1][k][d],u[l-2][k][d]);
               }
               numflux(UL, UR, 0.0, 1.0, flux);
               if(j==0)
               {
                  for(d=0; d<nvar; ++d)
                     res[j][i][d] -= dx * flux[d];
               }
               else if(j==nlocy)
               {
                  for(d=0; d<nvar; ++d)
                     res[j-1][i][d] += dx * flux[d];
               }
               else
               {
                  for(d=0; d<nvar; ++d)
                  {
                     res[j][i][d]   -= dx * flux[d];
                     res[j-1][i][d] += dx * flux[d];
                  }
               }
            }

         // Update solution
         for(j=jbeg; j<jbeg+nlocy; ++j)
            for(i=ibeg; i<ibeg+nlocx; ++i)
               for(d=0; d<nvar; ++d)
                  unew[j][i][d] = ark[rk]*uold[j-jbeg][i-ibeg][d]
                                 + (1.0-ark[rk])*(u[j][i][d] - lam * res[j-jbeg][i-ibeg][d]);

         ierr = DMDAVecRestoreArrayDOFRead(da, ul, &u); CHKERRQ(ierr);
         ierr = DMDAVecRestoreArrayDOF(da, ug, &unew); CHKERRQ(ierr);
      }

      t += dt; ++it;
      PetscPrintf(PETSC_COMM_WORLD,"it, t = %d, %f\n", it, t);
      if(it%si == 0 || PetscAbs(t-Tf) < 1.0e-13)
      {
         ierr = savesol(&c, t, da, ug); CHKERRQ(ierr);
      }
   }

   // Destroy everything before finishing
   ierr = VecDestroy(&ug); CHKERRQ(ierr);
   ierr = DMRestoreLocalVector(da, &ul); CHKERRQ(ierr);
   ierr = DMDestroy(&da); CHKERRQ(ierr);

   free(res); free(uold);

   ierr = PetscFinalize(); CHKERRQ(ierr);
}
