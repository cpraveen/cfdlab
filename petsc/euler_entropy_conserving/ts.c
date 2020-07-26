static char help[] = "Solves 2d Euler equations.\n\n";

#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>

#define min(a,b)  ( (a < b) ? a : b )
#define nvar  4

const PetscInt sw = 3; // stencil width, 3 on either side, for weno5
//const double xmin = -5.0, xmax = 5.0;
//const double ymin = -5.0, ymax = 5.0;
double xmin = -1.0, xmax = 1.0;
double ymin = -1.0, ymax = 1.0;
const double gas_gamma = 1.4;
const double gas_const = 1.0;
double dx, dy;

typedef enum { flux_central,flux_kepec2,flux_kepec4,flux_kep2,flux_mkep2,
               flux_mkep4,flux_kg2 } FluxScheme;
const char *const FluxSchemes[] = {"central","kepec2","kepec4","kep2","mkep2",
                                   "mkep4","kg2",
                                   "FluxScheme", "flux_", NULL};

typedef enum { prob_vortex, prob_density } Problem;
const char *const Problems[] = {"vortex", "density", "Problem", "prob_", NULL};

typedef struct
{
   PetscReal dt, cfl, Tf;
   PetscInt  max_steps, si;
   FluxScheme flux_scheme;
} AppCtx;

extern PetscErrorCode RHSFunction(TS,PetscReal,Vec,Vec,void*);
extern PetscErrorCode Monitor(TS,PetscInt,PetscReal,Vec,void*);

// Isentropic vortex
void initcond_vortex(const double x, const double y, double *Prim)
{
   const double M = 0.5;
   const double alpha = 45.0;
   const double beta = 5.0;
   const double r2 = x*x + y*y;
   Prim[0] =  pow(1.0 - (gas_gamma-1.0)*(beta*beta)/(8.0*gas_gamma*M_PI*M_PI)*exp(1-r2), (1.0/(gas_gamma-1.0)));
   Prim[1] =  M*cos(alpha*M_PI/180.0) - beta/(2.0*M_PI)*y*exp(0.5*(1.0-r2));
   Prim[2] =  M*sin(alpha*M_PI/180.0) + beta/(2.0*M_PI)*x*exp(0.5*(1.0-r2));
   Prim[3] =  pow(Prim[0],gas_gamma);
}

// Density profile
void initcond_density(const double x, const double y, double *Prim)
{
   Prim[0] = 1.0 + 0.98 * sin(2.0 * M_PI*(x + y));
   Prim[1] = 0.1;
   Prim[2] = 0.2;
   Prim[3] = 20.0;
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
   const double a = sqrt(gas_gamma*Prim[3]/Prim[0]);
   const double sx = fabs(Prim[1]) + a;
   const double sy = fabs(Prim[2]) + a;
   return 1.0/(sx/dx + sy/dy);
}

// Simple average flux
void avgflux(const double *Ul, const double *Ur, 
             const double nx, const double ny, 
             double *flux)
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

double logavg(double a, double b)
{
   double xi = b/a;
   double f = (xi - 1.0) / (xi + 1.0);
   double u = pow(f,2);

   double FF;
   if (u < 1.0e-2)
   {
      double u2 = pow(u,2);
      double u3 = u2 * u;
      FF = 1.0 + u/3.0 + u2/5.0 + u3/7.0;
   }
   else
      FF = 0.5 * log(xi) / f;

   return 0.5*(a+b)/FF;
}

// KEPEC flux
void numflux2(const double *Ul, const double *Ur,
              const double nx, const double ny,
              double *flux)
{
   double ql[nvar], qr[nvar];
   con2prim(Ul, ql);
   con2prim(Ur, qr);

   double ql2  = pow(ql[1],2) + pow(ql[2],2);
   double bl   = 0.5 * ql[0] / ql[3];

   double qr2  = pow(qr[1],2) + pow(qr[2],2);
   double br   = 0.5 * qr[0] / qr[3];

   double logr = logavg(ql[0], qr[0]);
   double logb = logavg(bl, br);
   double u    = 0.5*(ql[1] + qr[1]);
   double v    = 0.5*(ql[2] + qr[2]);
   double q2   = 0.5*(ql2 + qr2);

   double ra   = 0.5*(ql[0] + qr[0]);
   double ba   = 0.5*(bl + br);
   double p    = 0.5*ra/ba;

   // Rotated velocity
   double un   = u*nx + v*ny;

   // Centered flux
   flux[0] = logr*un;
   flux[1] = p*nx + u*flux[0];
   flux[2] = p*ny + v*flux[0];
   flux[3] = 0.5 * (1.0/((gas_gamma-1.0)*logb) - q2) * flux[0]
             + u*flux[1] + v*flux[2];
}

// 4th order KEPEC flux
void numflux4(const double *Ull, const double *Ul, 
              const double *Ur, const double *Urr,
              const double nx, const double ny,
              double *flux)
{
   double flux1[nvar], flux2[nvar], flux3[nvar];
   numflux2(Ul,Ur,nx,ny,flux1);
   numflux2(Ull,Ur,nx,ny,flux2);
   numflux2(Ul,Urr,nx,ny,flux3);
   for(int i=0; i<nvar; ++i)
      flux[i] = (4.0/3.0)*flux1[i] - (1.0/6.0)*flux2[i] - (1.0/6.0)*flux3[i];
}

// KEP flux of Jameson
void kepflux2(const double *Ul, const double *Ur,
               const double nx, const double ny,
               double *flux)
{
   double ql[nvar], qr[nvar];
   con2prim(Ul, ql);
   con2prim(Ur, qr);

   double ql2 = pow(ql[1], 2) + pow(ql[2], 2);
   double qr2 = pow(qr[1], 2) + pow(qr[2], 2);

   double al2 = gas_gamma * ql[3] / ql[0];
   double ar2 = gas_gamma * qr[3] / qr[0];

   double Hl = al2 / (gas_gamma-1.0) + 0.5 * ql2;
   double Hr = ar2 / (gas_gamma-1.0) + 0.5 * qr2;

   double r = 0.5 * (ql[0] + qr[0]);
   double u = 0.5 * (ql[1] + qr[1]);
   double v = 0.5 * (ql[2] + qr[2]);
   double p = 0.5 * (ql[3] + qr[3]);
   double H = 0.5 * (Hl + Hr);

   // Rotated velocity
   double un = u * nx + v * ny;

   // Centered flux
   flux[0] = r * un;
   flux[1] = p * nx + u * flux[0];
   flux[2] = p * ny + v * flux[0];
   flux[3] = r * H * un;
}
//------------------------------------------------------------------------------
// New KEP flux
void mkepflux2(const double *Ul, const double *Ur,
               const double nx, const double ny,
               double *flux)
{
   double ql[nvar], qr[nvar];
   con2prim(Ul, ql);
   con2prim(Ur, qr);

   double r = 0.5 * (ql[0] + qr[0]);
   double u = 0.5 * (ql[1] + qr[1]);
   double v = 0.5 * (ql[2] + qr[2]);
   double p = 0.5 * (ql[3] + qr[3]);
   double E = 0.5 * (Ul[3] + Ur[3]);

   // Rotated velocity
   double un = u * nx + v * ny;

   // Centered flux
   flux[0] = r * un;
   flux[1] = p * nx + u * flux[0];
   flux[2] = p * ny + v * flux[0];
   flux[3] = (E + p) * un;
}
// 4th order modified KEP flux
void mkepflux4(const double *Ull, const double *Ul,
               const double *Ur, const double *Urr,
               const double nx, const double ny,
               double *flux)
{
   double flux1[nvar], flux2[nvar], flux3[nvar];
   mkepflux2(Ul, Ur, nx, ny, flux1);
   mkepflux2(Ull, Ur, nx, ny, flux2);
   mkepflux2(Ul, Urr, nx, ny, flux3);
   for (int i = 0; i < nvar; ++i)
      flux[i] = (4.0 / 3.0) * flux1[i] - (1.0 / 6.0) * flux2[i] - (1.0 / 6.0) * flux3[i];
}
// Kennedy and Gruber flux
void kgflux2(const double *Ul, const double *Ur,
             const double nx, const double ny,
             double *flux)
{
   double ql[nvar], qr[nvar];
   con2prim(Ul, ql);
   con2prim(Ur, qr);

   double el = Ul[3]/Ul[0];
   double er = Ur[3]/Ur[0];

   double r = 0.5 * (ql[0] + qr[0]);
   double u = 0.5 * (ql[1] + qr[1]);
   double v = 0.5 * (ql[2] + qr[2]);
   double p = 0.5 * (ql[3] + qr[3]);
   double e = 0.5 * (el + er);

   // Rotated velocity
   double un = u * nx + v * ny;

   // Centered flux
   flux[0] = r * un;
   flux[1] = p * nx + u * flux[0];
   flux[2] = p * ny + v * flux[0];
   flux[3] = (r * e + p) * un;
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
   PetscReal   dtglobal, dtlocal = 1.0e20;
   PetscScalar ***u;
   Problem     problem;

   ierr = PetscInitialize(&argc, &argv, (char*)0, help); CHKERRQ(ierr);

   ctx.Tf  = 10.0;
   ctx.dt  = -1.0;
   ctx.cfl = -1.0;
   ctx.max_steps = 1000000;
   ctx.si = 100;

   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
   MPI_Comm_size(PETSC_COMM_WORLD, &size);
   
   // Get some command line options
   ierr = PetscOptionsGetReal(NULL,NULL,"-Tf",&ctx.Tf,NULL); CHKERRQ(ierr);
   ierr = PetscOptionsGetReal(NULL,NULL,"-dt",&ctx.dt,NULL); CHKERRQ(ierr);
   ierr = PetscOptionsGetReal(NULL,NULL,"-cfl",&ctx.cfl,NULL); CHKERRQ(ierr);
   ierr = PetscOptionsGetInt(NULL,NULL,"-si",&ctx.si,NULL); CHKERRQ(ierr);
   ierr = PetscOptionsGetEnum(NULL,NULL,"-flux",FluxSchemes,(PetscEnum *)&ctx.flux_scheme, NULL);
   ierr = PetscOptionsGetEnum(NULL,NULL,"-problem",Problems,(PetscEnum *)&problem, NULL);

   if(problem == prob_vortex)
   {
      xmin = -5.0;
      xmax = +5.0;
      ymin = -5.0;
      ymax = +5.0;
   }
   else if(problem == prob_density)
   {
      xmin = -1.0;
      xmax = +1.0;
      ymin = -1.0;
      ymax = +1.0;
   }
   else
   {
      PetscPrintf(PETSC_COMM_WORLD,"Unknown problem\n");
      exit(0);
   }

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
         if(problem == prob_vortex)
            initcond_vortex(x, y, prim);
         else if(problem == prob_density)
            initcond_density(x, y, prim);
         prim2con(prim, u[j][i]);
         dtlocal = min(dtlocal, dt_local(u[j][i]));
      }
   ierr = DMDAVecRestoreArrayDOF(da, ug, &u); CHKERRQ(ierr);
   MPI_Allreduce(&dtlocal, &dtglobal, 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);
   if(ctx.cfl > 0)
   {
      ctx.dt = ctx.cfl * dtglobal;
      PetscPrintf(PETSC_COMM_WORLD,"Using dt based on specified cfl = %f\n",ctx.cfl);
   }
   else if(ctx.dt > 0)
   {
      PetscPrintf(PETSC_COMM_WORLD,"Global dt = %e\n", dtglobal);
      PetscPrintf(PETSC_COMM_WORLD,"Using specified dt\n");
   }
   else
   {
      PetscPrintf(PETSC_COMM_WORLD,"Specify atleast dt or cfl\n");
      return(0);
   }
   PetscPrintf(PETSC_COMM_WORLD,"Initial time step = %e\n", ctx.dt);

   // Save initial condition to file
   ierr = savesol(0.0, da, ug); CHKERRQ(ierr);

   ierr = TSCreate(PETSC_COMM_WORLD,&ts); CHKERRQ(ierr);
   ierr = TSSetDM(ts,da); CHKERRQ(ierr);
   ierr = TSSetProblemType(ts,TS_NONLINEAR); CHKERRQ(ierr);
   ierr = TSSetRHSFunction(ts,NULL,RHSFunction,&ctx); CHKERRQ(ierr);
   ierr = TSSetTimeStep(ts,ctx.dt);
   ierr = TSSetType(ts,TSSSP); CHKERRQ(ierr);
   ierr = TSSetMaxSteps(ts,ctx.max_steps); CHKERRQ(ierr);
   ierr = TSSetMaxTime(ts,ctx.Tf); CHKERRQ(ierr);
   ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP); CHKERRQ(ierr);
   ierr = TSSetSolution(ts,ug); CHKERRQ(ierr);
   ierr = TSMonitorSet(ts,Monitor,&ctx,NULL); CHKERRQ(ierr);
   ierr = TSSetFromOptions(ts); CHKERRQ(ierr);
   ierr = TSSetUp(ts); CHKERRQ(ierr);

   ierr = TSSolve(ts,ug); CHKERRQ(ierr);

   // Destroy everything before finishing
   ierr = VecDestroy(&ug); CHKERRQ(ierr);
   ierr = DMDestroy(&da); CHKERRQ(ierr);
   ierr = TSDestroy(&ts); CHKERRQ(ierr);

   ierr = PetscFinalize(); CHKERRQ(ierr);
}

// The rhs function in du/dt = R(t,u)
PetscErrorCode RHSFunction(TS ts,PetscReal time,Vec U,Vec R,void* ptr)
{
   AppCtx*        ctx = (AppCtx*) ptr;
   DM             da;
   Vec            localU;
   PetscScalar    ***u;
   PetscScalar    ***res;
   PetscInt       i, j, ibeg, jbeg, nlocx, nlocy, d;
   PetscReal      flux[nvar], lam;
   PetscErrorCode ierr;

   ierr = TSGetDM(ts, &da); CHKERRQ(ierr);
   ierr = DMGetLocalVector(da,&localU); CHKERRQ(ierr);
   ierr = DMGlobalToLocalBegin(da, U, INSERT_VALUES, localU); CHKERRQ(ierr);
   ierr = DMGlobalToLocalEnd(da, U, INSERT_VALUES, localU); CHKERRQ(ierr);
   ierr = DMDAVecGetArrayDOFRead(da, localU, &u); CHKERRQ(ierr);
   ierr = DMDAVecGetArrayDOF(da, R, &res); CHKERRQ(ierr);

   ierr = DMDAGetCorners(da, &ibeg, &jbeg, 0, &nlocx, &nlocy, 0); CHKERRQ(ierr);

   // ---Begin res computation---
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
         if(ctx->flux_scheme == flux_central)
            avgflux(u[j][i-1], u[j][i], 1.0, 0.0, flux);
         else if(ctx->flux_scheme == flux_kepec2)
            numflux2(u[j][i-1], u[j][i], 1.0, 0.0, flux);
         else if(ctx->flux_scheme == flux_kepec4)
            numflux4(u[j][i-2], u[j][i-1], u[j][i], u[j][i+1], 1.0, 0.0, flux);
         else if(ctx->flux_scheme == flux_kep2)
            kepflux2(u[j][i-1], u[j][i], 1.0, 0.0, flux);
         else if(ctx->flux_scheme == flux_mkep2)
            mkepflux2(u[j][i-1], u[j][i], 1.0, 0.0, flux);
         else if(ctx->flux_scheme == flux_mkep4)
            mkepflux4(u[j][i-2], u[j][i-1], u[j][i], u[j][i+1], 1.0, 0.0, flux);
         else if(ctx->flux_scheme == flux_kg2)
            kgflux2(u[j][i-1], u[j][i], 1.0, 0.0, flux);
         else
         {
            PetscPrintf(PETSC_COMM_WORLD,"Unknown flux !!!\n");
            exit(0);
         }
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
         if(ctx->flux_scheme == flux_central)
            avgflux(u[j-1][i], u[j][i], 0.0, 1.0, flux);
         else if(ctx->flux_scheme == flux_kepec2)
            numflux2(u[j-1][i], u[j][i], 0.0, 1.0, flux);
         else if(ctx->flux_scheme == flux_kepec4)
            numflux4(u[j-2][i], u[j-1][i], u[j][i], u[j+1][i], 0.0, 1.0, flux);
         else if(ctx->flux_scheme == flux_kep2)
            kepflux2(u[j-1][i], u[j][i], 0.0, 1.0, flux);
         else if(ctx->flux_scheme == flux_mkep2)
            mkepflux2(u[j-1][i], u[j][i], 0.0, 1.0, flux);
         else if(ctx->flux_scheme == flux_mkep4)
            mkepflux4(u[j-2][i], u[j-1][i], u[j][i], u[j+1][i], 0.0, 1.0, flux);
         else if(ctx->flux_scheme == flux_kg2)
            kgflux2(u[j-1][i], u[j][i], 0.0, 1.0, flux);
         else
         {
            PetscPrintf(PETSC_COMM_WORLD,"Unknown flux !!!");
            exit(0);
         }
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
   // ---End res computation---

   ierr = DMDAVecRestoreArrayDOFRead(da, localU, &u); CHKERRQ(ierr);
   ierr = DMDAVecRestoreArrayDOF(da, R, &res); CHKERRQ(ierr);
   ierr = DMRestoreLocalVector(da,&localU); CHKERRQ(ierr);

   PetscFunctionReturn(0);
}

// This function is called after every time step.
PetscErrorCode Monitor(TS ts,PetscInt step,PetscReal time,Vec U,void *ptr)
{
   AppCtx*        ctx = (AppCtx*) ptr;
   DM             da;
   PetscInt       i, j, ibeg, jbeg, nlocx, nlocy;
   PetscReal      dtlocal, dtglobal;
   PetscScalar    ***u;
   PetscErrorCode ierr;

   if (step < 0) return(0); /* step of -1 indicates an interpolated solution */
   PetscPrintf(PETSC_COMM_WORLD,"iter, t = %d, %e\n", step, time);

   ierr = TSGetDM(ts, &da); CHKERRQ(ierr);

   if(step > 0 && (step%ctx->si == 0 || PetscAbs(time-ctx->Tf) < 1.0e-13))
   {
      ierr = savesol(time, da, U); CHKERRQ(ierr);
   }

   // If final time reached, dont do anything else, return from function.
   if(PetscAbs(time-ctx->Tf) < 1.0e-13)
      PetscFunctionReturn(0);

   // Compute time step based on cfl
   if(ctx->cfl > 0)
   {
      ierr = DMDAGetCorners(da, &ibeg, &jbeg, 0, &nlocx, &nlocy, 0); CHKERRQ(ierr);
      ierr = DMDAVecGetArrayDOFRead(da, U, &u); CHKERRQ(ierr);

      dtlocal = 1.0e20;
      for(j=jbeg; j<jbeg+nlocy; ++j)
         for(i=ibeg; i<ibeg+nlocx; ++i)
         {
            dtlocal = min(dtlocal, dt_local(u[j][i]));
         }
      ierr = DMDAVecRestoreArrayDOFRead(da, U, &u); CHKERRQ(ierr);
      MPI_Allreduce(&dtlocal, &dtglobal, 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);
      dtglobal *= ctx->cfl;
      // Adjust dt to reach final time exactly
      if(time+dtglobal > ctx->Tf) dtglobal = ctx->Tf - time;
      ierr = TSSetTimeStep(ts, dtglobal); CHKERRQ(ierr);
   }

   PetscFunctionReturn(0);
}
