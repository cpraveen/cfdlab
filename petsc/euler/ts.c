static char help[] = "Solves 2d Euler equations.\n\n";

#include <petscsys.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscvec.h>
#include <petscts.h>

#define min(a,b)  ( (a < b) ? a : b )
#define nvar  4

const double ark[] = {0.0, 3.0/4.0, 1.0/3.0};
const PetscInt sw = 3; // stencil width
const double xmin = -5.0, xmax = 5.0;
const double ymin = -5.0, ymax = 5.0;
const double gas_gamma = 1.4;
const double gas_const = 1.0;
double dx, dy;

typedef struct
{
   PetscReal dt, cfl, Tf;
   PetscInt  max_steps, si;
} AppCtx;

extern PetscErrorCode RHSFunction(TS,PetscReal,Vec,Vec,void*);
extern PetscErrorCode Monitor(TS,PetscInt,PetscReal,Vec,void*);

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

//------------------------------------------------------------------------------
// Weno reconstruction
// Return left value for face between u0, up1
//------------------------------------------------------------------------------
double weno5(const double um2, const double um1, const double u0, const double up1, const double up2)
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

// Conserved to primitive variables
void con2prim(const double *Con, double *Prim)
{
  Prim[0] = Con[0];
  Prim[1] = Con[1]/Con[0];
  Prim[2] = Con[2]/Con[0];
  Prim[3] = (Con[3] - 0.5*Prim[0]*(pow(Prim[1],2) + pow(Prim[2],2)))*(gas_gamma-1.0);
}

// Primitive to conserved variables
void prim2con(const double *Prim, double *Con)
{
  Con[0] = Prim[0];
  Con[1] = Prim[0]*Prim[1];
  Con[2] = Prim[0]*Prim[2];
  Con[3] = 0.5*Prim[0]*(pow(Prim[1],2) + pow(Prim[2],2)) + Prim[3]/(gas_gamma-1.0);
}

// Compute maximum eigenvalue in direction (nx,ny)
double maxeigval(const double *Con, const double nx, const double ny)
{
   double Prim[nvar];
   con2prim(Con, Prim);
   const double u = fabs(Prim[1]*nx + Prim[2]*ny);
   const double a = sqrt(gas_gamma*Prim[3]/Prim[0]);
   return u + a;
}

// Compute local timestep
double dt_local(const double *Con)
{
   double Prim[nvar];
   con2prim(Con, Prim);
   const double u = sqrt(pow(Prim[1],2) + pow(Prim[2],2));
   const double a = sqrt(gas_gamma*Prim[3]/Prim[0]);
   const double eigen = u + a;
   return min(dx,dy)/eigen;
}

// Simple average flux
void avg_flux(const double *Ul, const double *Ur, const double nx, const double ny, double *flux)
{
   double Pl[nvar], Pr[nvar];
   con2prim(Ul, Pl);
   con2prim(Ur, Pr);

   double fluxl[nvar], fluxr[nvar];

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
void numflux(const double *Ul, const double *Ur, const double nx, const double ny, double *flux)
{
   double Ua[nvar];
   for(int i=0; i<nvar; ++i) Ua[i] = 0.5*(Ul[i] + Ur[i]);
   double lam = maxeigval( Ua, nx, ny );
   avg_flux(Ul,Ur,nx,ny,flux);
   for(int i=0; i<nvar; ++i) flux[i] -= 0.5*lam*(Ur[i] - Ul[i]);
}

//------------------------------------------------------------------------------
PetscErrorCode savesol(double t, DM da, Vec ug)
{
   PetscErrorCode ierr;
   char           filename[32] = "sol";
   PetscMPIInt    rank;
   PetscInt       i, j, nx, ny, ibeg, jbeg, nlocx, nlocy;
   FILE           *fp;
   Vec            ul;
   PetscScalar    ***u;
   static int     c = 0;

   ierr = DMGetLocalVector(da, &ul); CHKERRQ(ierr);
   ierr = DMGlobalToLocalBegin(da, ug, INSERT_VALUES, ul); CHKERRQ(ierr);
   ierr = DMGlobalToLocalEnd(da, ug, INSERT_VALUES, ul); CHKERRQ(ierr);
   ierr = DMDAVecGetArrayDOFRead(da, ul, &u); CHKERRQ(ierr);
   ierr = DMDAGetInfo(da,0,&nx,&ny,0,0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);
   ierr = DMDAGetCorners(da, &ibeg, &jbeg, 0, &nlocx, &nlocy, 0); CHKERRQ(ierr);

   int iend = PetscMin(ibeg+nlocx+1, nx);
   int jend = PetscMin(jbeg+nlocy+1, ny);

   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
   sprintf(filename, "sol-%03d-%03d.plt", c, rank);
   fp = fopen(filename,"w");
   fprintf(fp, "TITLE = \"u_t + u_x + u_y = 0\"\n");
   fprintf(fp, "VARIABLES = x, y, rho, u, v, p\n");
   fprintf(fp, "ZONE STRANDID=1, SOLUTIONTIME=%e, I=%d, J=%d, DATAPACKING=POINT\n", t, iend-ibeg, jend-jbeg);
   for(j=jbeg; j<jend; ++j)
      for(i=ibeg; i<iend; ++i)
   {
      PetscReal x = xmin + i*dx + 0.5*dx;
      PetscReal y = ymin + j*dy + 0.5*dy;
      double prim[nvar];
      con2prim(u[j][i], prim);
      fprintf(fp, "%e %e %e %e %e %e\n", x, y, prim[0], prim[1], prim[2], prim[3]);
   }
   fclose(fp);

   ierr = DMDAVecRestoreArrayDOFRead(da, ul, &u); CHKERRQ(ierr);
   ierr = DMRestoreLocalVector(da, &ul); CHKERRQ(ierr);

   ++c;
   return(0);
}
//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
   // some parameters that can overwritten from command line
   PetscInt  nx  = 50, ny=50; // use -da_grid_x, -da_grid_y to override these
   
   PetscErrorCode ierr;
   AppCtx      ctx;
   TS          ts;
   DM          da;
   Vec         ug;
   PetscInt    i, j, ibeg, jbeg, nlocx, nlocy;
   PetscMPIInt rank, size;
   PetscScalar ***u;

   ierr = PetscInitialize(&argc, &argv, (char*)0, help); CHKERRQ(ierr);

   ctx.Tf  = 10.0;
   ctx.cfl = 0.4;
   ctx.max_steps = 1000000;
   ctx.si = 100;

   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
   MPI_Comm_size(PETSC_COMM_WORLD, &size);
   
   // Get some command line options
   ierr = PetscOptionsGetReal(NULL,NULL,"-Tf",&ctx.Tf,NULL); CHKERRQ(ierr);
   ierr = PetscOptionsGetReal(NULL,NULL,"-cfl",&ctx.cfl,NULL); CHKERRQ(ierr);
   ierr = PetscOptionsGetInt(NULL,NULL,"-si",&ctx.si,NULL); CHKERRQ(ierr);

   ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC,
                       DMDA_STENCIL_BOX, -nx, -ny, PETSC_DECIDE, PETSC_DECIDE, nvar,
                       sw, NULL, NULL, &da); CHKERRQ(ierr);
   ierr = DMDAGetInfo(da,0,&nx,&ny,0,0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);
   dx = (xmax - xmin) / (PetscReal)(nx);
   dy = (ymax - ymin) / (PetscReal)(ny);
   PetscPrintf(PETSC_COMM_WORLD,"nx = %d, dx = %e\n", nx, dx);
   PetscPrintf(PETSC_COMM_WORLD,"ny = %d, dy = %e\n", ny, dy);

   ierr = DMCreateGlobalVector(da, &ug); CHKERRQ(ierr);
   ierr = PetscObjectSetName((PetscObject) ug, "Solution"); CHKERRQ(ierr);

   ierr = DMDAGetCorners(da, &ibeg, &jbeg, 0, &nlocx, &nlocy, 0); CHKERRQ(ierr);
   ierr = DMDAVecGetArrayDOF(da, ug, &u); CHKERRQ(ierr);
   double dtlocal = 1.0e20;
   for(j=jbeg; j<jbeg+nlocy; ++j)
      for(i=ibeg; i<ibeg+nlocx; ++i)
      {
         PetscReal x = xmin + i*dx + 0.5*dx;
         PetscReal y = ymin + j*dy + 0.5*dy;
         double prim[nvar];
         initcond(x, y, prim);
         prim2con(prim, u[j][i]);
         dtlocal = min(dtlocal, dt_local(u[j][i]));
      }
   ierr = DMDAVecRestoreArrayDOF(da, ug, &u); CHKERRQ(ierr);
   MPI_Allreduce(&dtlocal, &ctx.dt, 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);
   ctx.dt *= ctx.cfl;
   PetscPrintf(PETSC_COMM_WORLD,"Initial time step = %e\n", ctx.dt);

   // Save initial condition to file
   ierr = savesol(0.0, da, ug); CHKERRQ(ierr);

   ierr = TSCreate(PETSC_COMM_WORLD,&ts); CHKERRQ(ierr);
   ierr = TSSetDM(ts,da); CHKERRQ(ierr);
   ierr = TSSetProblemType(ts,TS_NONLINEAR); CHKERRQ(ierr);
   ierr = TSSetRHSFunction(ts,NULL,RHSFunction,&ctx); CHKERRQ(ierr);
   ierr = TSSetInitialTimeStep(ts,0.0,ctx.dt);
   ierr = TSSetType(ts,TSSSP); CHKERRQ(ierr);
   ierr = TSSetDuration(ts,ctx.max_steps,ctx.Tf); CHKERRQ(ierr);
   ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP); CHKERRQ(ierr);
   ierr = TSSetSolution(ts,ug); CHKERRQ(ierr);
   ierr = TSMonitorSet(ts,Monitor,&ctx,NULL); CHKERRQ(ierr);
   ierr = TSSetFromOptions(ts); CHKERRQ(ierr);

   ierr = TSSolve(ts,ug); CHKERRQ(ierr);

   // Destroy everything before finishing
   ierr = VecDestroy(&ug); CHKERRQ(ierr);
   ierr = DMDestroy(&da); CHKERRQ(ierr);
   ierr = TSDestroy(&ts); CHKERRQ(ierr);

   ierr = PetscFinalize(); CHKERRQ(ierr);
}

PetscErrorCode RHSFunction(TS ts,PetscReal time,Vec U,Vec R,void* ptr)
{
   AppCtx*        ctx = (AppCtx*) ptr;
   DM             da;
   Vec            localU;
   PetscScalar    ***u;
   PetscScalar    ***res;
   PetscInt       i, j, ibeg, jbeg, nlocx, nlocy, d;
   PetscInt       il, jl, nl, ml;
   PetscReal      UL[nvar], UR[nvar], flux[nvar], lam;
   PetscErrorCode ierr;

   ierr = TSGetDM(ts, &da); CHKERRQ(ierr);
   ierr = DMGetLocalVector(da,&localU); CHKERRQ(ierr);
   ierr = DMGlobalToLocalBegin(da, U, INSERT_VALUES, localU); CHKERRQ(ierr);
   ierr = DMGlobalToLocalEnd(da, U, INSERT_VALUES, localU); CHKERRQ(ierr);
   ierr = DMDAVecGetArrayDOFRead(da, localU, &u); CHKERRQ(ierr);
   ierr = DMDAVecGetArrayDOF(da, R, &res); CHKERRQ(ierr);

   ierr = DMDAGetCorners(da, &ibeg, &jbeg, 0, &nlocx, &nlocy, 0); CHKERRQ(ierr);
   ierr = DMDAGetGhostCorners(da,&il,&jl,0,&nl,&ml,0); CHKERRQ(ierr);

   // Set residual 0
   for(j=jbeg; j<jbeg+nlocy; ++j)
      for(i=ibeg; i<ibeg+nlocx; ++i)
         for(d=0; d<nvar; ++d)
            res[j][i][d] = 0;

   // x fluxes
   for(i=ibeg; i<ibeg+nlocx+1; ++i)
      for(j=jbeg; j<jbeg+nlocy; ++j)
      {
         // face between i-1, i
         for(d=0; d<nvar; ++d)
         {
            UL[d] = weno5(u[j][i-3][d],u[j][i-2][d],u[j][i-1][d],u[j][i][d],u[j][i+1][d]);
            UR[d] = weno5(u[j][i+2][d],u[j][i+1][d],u[j][i][d],u[j][i-1][d],u[j][i-2][d]);
         }
         numflux(UL, UR, 1.0, 0.0, flux);
         if(i==ibeg)
         {
            for(d=0; d<nvar; ++d)
               res[j][i][d] -= dy * flux[d];
         }
         else if(i==ibeg+nlocx)
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
   for(j=jbeg; j<jbeg+nlocy+1; ++j)
      for(i=ibeg; i<ibeg+nlocx; ++i)
      {
         // face between j-1, j
         for(d=0; d<nvar; ++d)
         {
            UL[d] = weno5(u[j-3][i][d],u[j-2][i][d],u[j-1][i][d],u[j][i][d],u[j+1][i][d]);
            UR[d] = weno5(u[j+2][i][d],u[j+1][i][d],u[j][i][d],u[j-1][i][d],u[j-2][i][d]);
         }
         numflux(UL, UR, 0.0, 1.0, flux);
         if(j==jbeg)
         {
            for(d=0; d<nvar; ++d)
               res[j][i][d] -= dx * flux[d];
         }
         else if(j==jbeg+nlocy)
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

   lam = 1.0/(dx*dy);
   for(j=jbeg; j<jbeg+nlocy; ++j)
      for(i=ibeg; i<ibeg+nlocx; ++i)
         for(d=0; d<nvar; ++d)
            res[j][i][d] *= -lam;

   ierr = DMDAVecRestoreArrayDOFRead(da, localU, &u); CHKERRQ(ierr);
   ierr = DMDAVecRestoreArrayDOF(da, R, &res); CHKERRQ(ierr);
   ierr = DMRestoreLocalVector(da,&localU); CHKERRQ(ierr);

   PetscFunctionReturn(0);
}

PetscErrorCode Monitor(TS ts,PetscInt step,PetscReal time,Vec u,void *ptr)
{
   AppCtx*        ctx = (AppCtx*) ptr;
   DM             da;
   PetscErrorCode ierr;

   if (step < 0) return(0); /* step of -1 indicates an interpolated solution */

   if(step%ctx->si == 0 || PetscAbs(time-ctx->Tf) < 1.0e-13)
   {
      ierr = TSGetDM(ts, &da); CHKERRQ(ierr);
      ierr = savesol(time, da, u); CHKERRQ(ierr);
   }

   return(0);
}
