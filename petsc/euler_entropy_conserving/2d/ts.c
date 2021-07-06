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

typedef enum { flux_central,flux_kepec,flux_kep,flux_mkep,
               flux_kg,flux_ducros,flux_mkepec,flux_keep } FluxScheme;
const char *const FluxSchemes[] = {"central","kepec","kep","mkep",
                                   "kg","ducros","mkepec","keep",
                                   "FluxScheme", "flux_", NULL};

typedef enum { prob_vortex, prob_density } Problem;
const char *const Problems[] = {"vortex", "density", "Problem", "prob_", NULL};

typedef struct
{
   PetscReal dt, cfl, Tf;
   PetscInt  order, max_steps, si;
   FluxScheme flux_scheme;
} AppCtx;

extern PetscErrorCode RHSFunction(TS,PetscReal,Vec,Vec,void*);
extern PetscErrorCode Monitor(TS,PetscInt,PetscReal,Vec,void*);

// Isentropic vortex
void initcond_vortex(const double x, const double y, double *Prim)
{
   const double M = 0.5;
   const double alpha = M_PI*(45.0/180.0);
   const double beta = 5.0;
   const double r2 = x*x + y*y;
   const double T =  1.0 - (gas_gamma-1.0)*(beta*beta)/(8.0*gas_gamma*M_PI*M_PI)*exp(1-r2);
   Prim[0] =  pow(T, 1.0/(gas_gamma-1.0));
   Prim[1] =  M*cos(alpha) - beta/(2.0*M_PI)*y*exp(0.5*(1.0-r2));
   Prim[2] =  M*sin(alpha) + beta/(2.0*M_PI)*x*exp(0.5*(1.0-r2));
   Prim[3] =  pow(Prim[0],gas_gamma);
}

// Density profile
const double vx0 = 0.1;
const double vy0 = 0.2;
const double p0 = 20.0;
void initcond_density(const double x, const double y, double *Prim)
{
   const double eps = 0.0; // perturb velocity
   Prim[0] = 1.0 + 0.98 * sin(2.0 * M_PI * (x + y));
   Prim[1] = vx0 + eps * (sin(2 * M_PI * x) + sin(2 * M_PI * y));
   Prim[2] = vy0 + eps * (cos(2 * M_PI * x) + cos(2 * M_PI * y));
   Prim[3] = p0;
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
void kepecflux(const double *Ul, const double *Ur,
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

// Modified KEPEC flux
void mkepecflux(const double *Ul, const double *Ur,
                const double nx, const double ny,
                double *flux)
{
   double ql[nvar], qr[nvar];
   con2prim(Ul, ql);
   con2prim(Ur, qr);

   double ql2 = pow(ql[1], 2) + pow(ql[2], 2);
   double bl = 0.5 * ql[0] / ql[3];

   double qr2 = pow(qr[1], 2) + pow(qr[2], 2);
   double br = 0.5 * qr[0] / qr[3];

   double u = 0.5 * (ql[1] + qr[1]);
   double v = 0.5 * (ql[2] + qr[2]);
   double q2 = 0.5 * (ql2 + qr2);

   double ra = 0.5 * (ql[0] + qr[0]);
   double ba = 0.5 * (bl + br);

   // Rotated velocity
   double un = u * nx + v * ny;

   // Relative pressure difference
   const double r1 = 0.005;
   const double r2 = 0.01;
   double dp = 2.0*fabs(qr[3] - ql[3])/(ql[3] + qr[3]);
   double p, logr, logb;
   if(dp < r1) // use arithmetic average
   {
      p = 0.5*(ql[3] + qr[3]);
      logb = ba;
      logr = ra;
   }
   else
   {
      logr = logavg(ql[0], qr[0]);
      logb = logavg(bl, br);
      p = 0.5 * ra / ba;
      if(dp < r2) // r1 <= dp <= r2: blend between average and KEPEC
      {
         double x = (dp - r1)/(r2 - r1);
         double theta = 0.5*(1.0 + cos(M_PI*(x-1.0)));
         double p1 = 0.5*(ql[3] + qr[3]);
         p = theta * p + (1.0 - theta) * p1;
         logr = theta * logr + (1.0 - theta)* ra;
         logb = theta * logb + (1.0 - theta)* ba;
      }
   }

   // Centered flux
   flux[0] = logr * un;
   flux[1] = p * nx + u * flux[0];
   flux[2] = p * ny + v * flux[0];
   flux[3] = 0.5 * (1.0 / ((gas_gamma - 1.0) * logb) - q2) * flux[0] + u * flux[1] + v * flux[2];
}

// KEP flux of Jameson
void kepflux(const double *Ul, const double *Ur,
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
void mkepflux(const double *Ul, const double *Ur,
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
   double u2l = ql[1] * ql[1] + ql[2] * ql[2];
   double u2r = qr[1] * qr[1] + qr[2] * qr[2];
   double u2  = 0.5 * (u2l + u2r);

   // Rotated velocity
   double un = u * nx + v * ny;

   // Centered flux
   flux[0] = r * un;
   flux[1] = p * nx + u * flux[0];
   flux[2] = p * ny + v * flux[0];
   flux[3] = gas_gamma * p * un / (gas_gamma-1.0) + 0.5 * r * u2 * un;
}
//------------------------------------------------------------------------------
// Kennedy and Gruber flux
void kgflux(const double *Ul, const double *Ur,
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
// Ducros flux
void ducrosflux(const double *Ul, const double *Ur,
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
   double ru= 0.5 * (Ul[1] + Ur[1]);
   double rv= 0.5 * (Ul[2] + Ur[2]);
   double E = 0.5 * (Ul[3] + Ur[3]);

   // Rotated velocity
   double un = u * nx + v * ny;

   // Centered flux
   flux[0] = r * un;
   flux[1] = p * nx + ru * un;
   flux[2] = p * ny + rv * un;
   flux[3] = (E + p) * un;
}
//------------------------------------------------------------------------------
// Shima et al. flux: KEEP-PE
void keepflux(const double *Ul, const double *Ur,
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
   double k = 0.5 * (ql[1] * qr[1] + ql[2] * qr[2]);

   // Rotated velocity
   double unl = ql[1] * nx + ql[2] * ny;
   double unr = qr[1] * nx + qr[2] * ny;
   double un = 0.5 * (unl + unr);
   double pv = 0.5 * (ql[3] * unr + qr[3] * unl);

   // Centered flux
   flux[0] = r * un;
   flux[1] = p * nx + u * flux[0];
   flux[2] = p * ny + v * flux[0];
   flux[3] = p * un / (gas_gamma - 1.0) + r * k * un + pv;
}
//------------------------------------------------------------------------------
void numflux(const FluxScheme flux_scheme,
             const int order,
             const double *Ull, const double *Ul,
             const double *Ur, const double *Urr,
             const double nx, const double ny,
             double *flux)
{
   if(order == 2)
   {
      if (flux_scheme == flux_central)
         avgflux(Ul, Ur, nx, ny, flux);
      else if (flux_scheme == flux_kep)
         kepflux(Ul, Ur, nx, ny, flux);
      else if (flux_scheme == flux_kepec)
         kepecflux(Ul, Ur, nx, ny, flux);
      else if (flux_scheme == flux_mkep)
         mkepflux(Ul, Ur, nx, ny, flux);
      else if (flux_scheme == flux_kg)
         kgflux(Ul, Ur, nx, ny, flux);
      else if (flux_scheme == flux_ducros)
         ducrosflux(Ul, Ur, nx, ny, flux);
      else if (flux_scheme == flux_mkepec)
         mkepecflux(Ul, Ur, nx, ny, flux);
      else if (flux_scheme == flux_keep)
         keepflux(Ul, Ur, nx, ny, flux);
      else
      {
         PetscPrintf(PETSC_COMM_WORLD,"numflux: flux is not implemented\n");
         exit(0);
      }
   }
   else if(order == 4)
   {
      double flux1[nvar], flux2[nvar], flux3[nvar];
      if (flux_scheme == flux_central)
      {
         avgflux(Ul, Ur, nx, ny, flux1);
         avgflux(Ull, Ur, nx, ny, flux2);
         avgflux(Ul, Urr, nx, ny, flux3);
      }
      else if (flux_scheme == flux_kep)
      {
         kepflux(Ul, Ur, nx, ny, flux1);
         kepflux(Ull, Ur, nx, ny, flux2);
         kepflux(Ul, Urr, nx, ny, flux3);
      }
      else if (flux_scheme == flux_kepec)
      {
         kepecflux(Ul, Ur, nx, ny, flux1);
         kepecflux(Ull, Ur, nx, ny, flux2);
         kepecflux(Ul, Urr, nx, ny, flux3);
      }
      else if (flux_scheme == flux_mkep)
      {
         mkepflux(Ul, Ur, nx, ny, flux1);
         mkepflux(Ull, Ur, nx, ny, flux2);
         mkepflux(Ul, Urr, nx, ny, flux3);
      }
      else if (flux_scheme == flux_kg)
      {
         kgflux(Ul, Ur, nx, ny, flux1);
         kgflux(Ull, Ur, nx, ny, flux2);
         kgflux(Ul, Urr, nx, ny, flux3);
      }
      else if (flux_scheme == flux_ducros)
      {
         ducrosflux(Ul, Ur, nx, ny, flux1);
         ducrosflux(Ull, Ur, nx, ny, flux2);
         ducrosflux(Ul, Urr, nx, ny, flux3);
      }
      else if (flux_scheme == flux_mkepec)
      {
         mkepecflux(Ul, Ur, nx, ny, flux1);
         mkepecflux(Ull, Ur, nx, ny, flux2);
         mkepecflux(Ul, Urr, nx, ny, flux3);
      }
      else if (flux_scheme == flux_keep)
      {
         keepflux(Ul, Ur, nx, ny, flux1);
         keepflux(Ull, Ur, nx, ny, flux2);
         keepflux(Ul, Urr, nx, ny, flux3);
      }
      else
      {
         PetscPrintf(PETSC_COMM_WORLD,"numflux: flux is not implemented\n");
         exit(0);
      }
      for (int i = 0; i < nvar; ++i)
         flux[i] = (4.0 / 3.0) * flux1[i] - (1.0 / 6.0) * flux2[i] - (1.0 / 6.0) * flux3[i];
   }
   else
   {
      PetscPrintf(PETSC_COMM_WORLD,"numflux: order is not implemented\n");
      exit(0);
   }
}

//------------------------------------------------------------------------------
PetscErrorCode compute_global(PetscReal time, DM da, Vec ug, PetscReal *global)
{
   PetscErrorCode ierr;
   PetscScalar    ***u;
   PetscInt       i, j, ibeg, jbeg, nlocx, nlocy;
   PetscReal      local[2];
   ierr = DMDAVecGetArrayDOFRead(da, ug, &u); CHKERRQ(ierr);
   ierr = DMDAGetCorners(da, &ibeg, &jbeg, 0, &nlocx, &nlocy, 0); CHKERRQ(ierr);
   for(i=0; i<2; ++i)
      local[i] = 0.0;
   for(j=jbeg; j<jbeg+nlocy; ++j)
      for(i=ibeg; i<ibeg+nlocx; ++i)
      {
         double prim[nvar];
         con2prim(u[j][i], prim);
         if(0) // set to 1 for density wave test
         {
            PetscReal x = xmin + i * dx + 0.5 * dx;
            PetscReal y = ymin + j * dy + 0.5 * dy;
            double prim0[nvar];
            initcond_density(x-vx0*time, y-vy0*time, prim0);
            double dp = prim[3] - p0;
            double dvx = prim[1] - vx0;
            double dvy = prim[2] - vy0;
            local[0] += 0.5*dp*dp/p0/gas_gamma + 0.5*prim0[0]*(dvx*dvx + dvy*dvy);
         }
         else
         {
            local[0] += 0.5 * prim[0] * (pow(prim[1],2) + pow(prim[2],2));
         }
         local[1] += -prim[0] * (log(prim[3]) - gas_gamma * log(prim[0]));
      }
   for(i=0; i<2; ++i)
      local[i] *= dx * dy;
   // sum over all procs
   MPI_Reduce(local, global, 2, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
   ierr = DMDAVecRestoreArrayDOFRead(da, ug, &u); CHKERRQ(ierr);
   PetscFunctionReturn(0);
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
   ierr = PetscOptionsGetEnum(NULL,NULL,"-flux",FluxSchemes,
                              (PetscEnum *)&ctx.flux_scheme, NULL); CHKERRQ(ierr);
   ierr = PetscOptionsGetInt(NULL,NULL,"-order",&ctx.order,NULL); CHKERRQ(ierr);
   ierr = PetscOptionsGetEnum(NULL,NULL,"-problem",Problems,
                              (PetscEnum *)&problem, NULL); CHKERRQ(ierr);

   if(problem == prob_vortex)
   {
      xmin = -10.0;
      xmax = +10.0;
      ymin = -10.0;
      ymax = +10.0;
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
         numflux(ctx->flux_scheme, ctx->order, u[j][i-2], u[j][i-1], u[j][i], u[j][i+1], 1.0, 0.0, flux);
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
         numflux(ctx->flux_scheme, ctx->order, u[j-2][i], u[j-1][i], u[j][i], u[j+1][i], 0.0, 1.0, flux);
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
   PetscReal      global[2], dtlocal, dtglobal;
   PetscScalar    ***u;
   PetscErrorCode ierr;

   if (step < 0) return(0); /* step of -1 indicates an interpolated solution */

   ierr = TSGetDM(ts, &da); CHKERRQ(ierr);
   ierr = compute_global(time, da, U, global); CHKERRQ(ierr);
   PetscPrintf(PETSC_COMM_WORLD,"it,t,ke,ent= %d %18.10e %18.10e %18.10e\n", step, time, global[0], global[1]);

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
