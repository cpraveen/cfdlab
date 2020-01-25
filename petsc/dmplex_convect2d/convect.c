static char help[] = "Linear advection using first order FVM\n";
#include <petsc.h>

PETSC_STATIC_INLINE
void Waxpy(PetscInt dim, PetscScalar a, const PetscScalar *x, const PetscScalar *y, PetscScalar *w)
{
   PetscInt d; for (d = 0; d < dim; ++d) w[d] = a*x[d] + y[d];
}

PETSC_STATIC_INLINE
PetscScalar Dot(PetscInt dim, const PetscScalar *x, const PetscScalar *y)
{
   PetscScalar sum = 0.0; PetscInt d; for (d = 0; d < dim; ++d) sum += x[d]*y[d]; return sum;
}

PETSC_STATIC_INLINE
PetscReal Norm(PetscInt dim, const PetscScalar *x)
{
   return PetscSqrtReal(PetscAbsScalar(Dot(dim,x,x)));
}

typedef struct
{
   PetscReal *X0,*V0;     // 0-height centroids and volumes
   PetscReal *X1,*V1,*N1; // 1-height centroids, volumes, and normals
} MeshInfo;

typedef struct
{
   PetscInt dim;
   MeshInfo info;
   char     outputBasename[PETSC_MAX_PATH_LEN];
   PetscInt vtkInterval;
} AppCtx;

void boundary_value(const PetscReal *X, PetscReal *u)
{
   *u = 0.0;
}

void advection_speed(const PetscReal *X, PetscReal *speed)
{
   PetscReal x = X[0];
   PetscReal y = X[1];
   speed[0] = -y;
   speed[1] =  x;
}

void initial_condition(const PetscReal *X,
                       PetscReal       *u)
{
   PetscReal x = X[0];
   PetscReal y = X[1];
   PetscReal r2 = PetscPowReal(x-0.5,2) + PetscPowReal(y,2);
   *u = PetscExpReal(-50.0*r2);
}

// Upwind flux
void numerical_flux(const PetscReal a,
                    const PetscReal *ul,
                    const PetscReal *ur,
                    PetscReal       *Flux)
{
   *Flux = (a > 0.0) ? a*(*ul) : a*(*ur);
}

#undef __FUNCT__
#define __FUNCT__ "MeshInfoCreate"
PetscErrorCode MeshInfoCreate(DM dm,MeshInfo *info)
{
   PetscFunctionBegin;
   PetscErrorCode ierr;
   PetscInt       dim,i,p,pStart,pEnd,nCells;
   PetscReal      dummy[3];
   ierr   = DMGetDimension(dm,&dim);CHKERRQ(ierr);

   // cell geometry info
   ierr   = DMPlexGetHeightStratum(dm,0,&pStart,&pEnd);CHKERRQ(ierr);
   nCells = pEnd-pStart;
   ierr   = PetscMalloc2(nCells*dim*sizeof(PetscReal),&(info->X0),
                         nCells    *sizeof(PetscReal),&(info->V0));CHKERRQ(ierr);
   for(p=pStart;p<pEnd;p++)
   {
      i    = p-pStart;
      ierr = DMPlexComputeCellGeometryFVM(dm,p,&(info->V0[i]),&(info->X0[i*dim]),dummy);CHKERRQ(ierr);
   }

   // face geometry info
   ierr   = DMPlexGetHeightStratum(dm,1,&pStart,&pEnd);CHKERRQ(ierr);
   nCells = pEnd-pStart;
   ierr   = PetscMalloc3(nCells*dim*sizeof(PetscReal),&(info->X1),
                         nCells*dim*sizeof(PetscReal),&(info->N1),
                         nCells*    sizeof(PetscReal),&(info->V1));CHKERRQ(ierr);
   for(p=pStart;p<pEnd;p++)
   {
      i    = p-pStart;
      ierr = DMPlexComputeCellGeometryFVM(dm,p,&(info->V1[i]),&(info->X1[i*dim]),&(info->N1[i*dim]));CHKERRQ(ierr);
   }

   PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MeshInfoDestroy"
PetscErrorCode MeshInfoDestroy(MeshInfo *info)
{
   PetscFunctionBegin;
   PetscErrorCode ierr;
   ierr = PetscFree2(info->X0,info->V0);CHKERRQ(ierr);
   ierr = PetscFree3(info->X1,info->V1,info->N1);CHKERRQ(ierr);
   PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SetIC"
PetscErrorCode SetIC(DM dm, void *ctx, Vec U)
{
   PetscFunctionBegin;
   AppCtx        *user  = (AppCtx *)ctx;
   MeshInfo      *info  = &(user->info);
   PetscErrorCode ierr;
   PetscInt dim = user->dim;

   PetscScalar   *Uarray;
   ierr = VecGetArray(U, &Uarray); CHKERRQ(ierr);

   PetscInt p, pStart, pEnd;
   ierr = DMPlexGetHeightStratum(dm,0,&pStart,&pEnd);CHKERRQ(ierr);
   for(p=pStart;p<pEnd;p++)
   {
      PetscReal *X = &(info->X0[(p-pStart)*dim]);
      PetscScalar uinit, *u;
      initial_condition(X, &uinit);
      ierr = DMPlexPointGlobalRef(dm,p-pStart,Uarray,&u);CHKERRQ(ierr);
      *u = uinit;
   }

   ierr = VecRestoreArray(U, &Uarray); CHKERRQ(ierr);
   PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "OutputVTK"
static PetscErrorCode OutputVTK(DM dm, const char *filename, PetscViewer *viewer)
{
   PetscFunctionBegin;
   PetscErrorCode ierr;

   ierr = PetscViewerCreate(PetscObjectComm((PetscObject)dm), viewer);CHKERRQ(ierr);
   ierr = PetscViewerSetType(*viewer, PETSCVIEWERVTK);CHKERRQ(ierr);
   ierr = PetscViewerFileSetName(*viewer, filename);CHKERRQ(ierr);
   PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Monitor"
PetscErrorCode Monitor(TS ts,PetscInt step,PetscReal time,Vec U,void *ctx)
{
   PetscFunctionBegin;
   AppCtx         *user  = (AppCtx *)ctx;
   DM             dm;
   PetscErrorCode ierr;
   static PetscInt counter = 0;

   if (step < 0) return(0); /* step of -1 indicates an interpolated solution */

   if(step % user->vtkInterval != 0) return(0);

   char filename[PETSC_MAX_PATH_LEN];
   ierr = PetscSNPrintf(filename,sizeof filename,"%s-%03D.vtk",user->outputBasename,counter); CHKERRQ(ierr);
   ierr = TSGetDM(ts, &dm); CHKERRQ(ierr);
   PetscViewer viewer;
   ierr = OutputVTK(dm, filename, &viewer); CHKERRQ(ierr);
   ierr = VecView(U, viewer); CHKERRQ(ierr);
   ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
   ++counter;
   PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeResidual"
PetscErrorCode ComputeResidual(TS ts,PetscReal t,Vec U,Vec R,void *ctx)
{
   PetscFunctionBegin;
   AppCtx        *user  = (AppCtx *)ctx;
   MeshInfo      *info  = &(user->info);
   PetscErrorCode ierr;
   DM             dm;
   ierr = TSGetDM(ts,&dm);CHKERRQ(ierr);

   // Global --> Local of the solution
   Vec          Ul;
   PetscScalar *Pl;
   ierr = DMGetLocalVector(dm,&Ul  );CHKERRQ(ierr);
   ierr = DMGlobalToLocalBegin(dm,U  ,INSERT_VALUES,Ul  );CHKERRQ(ierr);
   ierr = DMGlobalToLocalEnd  (dm,U  ,INSERT_VALUES,Ul  );CHKERRQ(ierr);
   ierr = VecGetArray(Ul  ,&Pl  );CHKERRQ(ierr);

   // Set rhs to zero. Is this really required ?
   ierr = VecSet(R,0.0);CHKERRQ(ierr);

   // Compute the residual
   PetscInt     p,pStart,pEnd,ss,dim = user->dim;
   PetscReal   *Rarray,*r,*ul,*ur,Flux,speed[dim],normal_speed;
   ierr = VecGetArray(R,&Rarray);CHKERRQ(ierr);

   // Loop over faces
   ierr = DMPlexGetHeightStratum(dm,1,&pStart,&pEnd);CHKERRQ(ierr);
   for(p=pStart;p<pEnd;p++)
   {
      // Get support of interface
      const PetscInt *supp;
      ierr = DMPlexGetSupportSize(dm,p,&ss  );CHKERRQ(ierr);
      ierr = DMPlexGetSupport    (dm,p,&supp);CHKERRQ(ierr);

      PetscReal *Xf = &(info->X1[(p-pStart)*dim]); // face mid-point
      PetscReal *Nf = &(info->N1[(p-pStart)*dim]); // face normal
      advection_speed(Xf, speed);
      normal_speed = Dot(dim, speed, Nf);

      if(ss==1) // Boundary face
      {
         PetscReal *Xl = &(info->X0[supp[0]*dim]);

         //ierr = DMPlexPointGlobalRef(dm,supp[0],Pl,&ul);CHKERRQ(ierr);
         ul = &Pl[supp[0]];
         boundary_value(Xf, ur);
         numerical_flux(normal_speed, ul, ur, &Flux);
         ierr = DMPlexPointGlobalRef(dm,supp[0],Rarray,&r);CHKERRQ(ierr);
         if(r) *r -= Flux * info->V1[p-pStart] / info->V0[supp[0]];
      }
      else // Interior face
      {
         PetscReal *Xl = &(info->X0[supp[0]*dim]);
         PetscReal *Xr = &(info->X0[supp[1]*dim]);

         //ierr = DMPlexPointGlobalRef(dm,supp[0],Pl,&ul);CHKERRQ(ierr);
         //ierr = DMPlexPointGlobalRef(dm,supp[1],Pl,&ur);CHKERRQ(ierr);
         ul = &Pl[supp[0]];
         ur = &Pl[supp[1]];
         numerical_flux(normal_speed, ul, ur, &Flux);
         ierr = DMPlexPointGlobalRef(dm,supp[0],Rarray,&r);CHKERRQ(ierr);
         if(r) *r -= Flux * info->V1[p-pStart] / info->V0[supp[0]];
         ierr = DMPlexPointGlobalRef(dm,supp[1],Rarray,&r);CHKERRQ(ierr);
         if(r) *r += Flux * info->V1[p-pStart] / info->V0[supp[1]];
      }
   }

   // Cleanup
   ierr = VecRestoreArray(Ul  ,&Pl    );CHKERRQ(ierr);
   ierr = VecRestoreArray(R   ,&Rarray);CHKERRQ(ierr);
   ierr = DMRestoreLocalVector(dm,&Ul  );CHKERRQ(ierr);
   PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
   // Initialize
   MPI_Comm          comm;
   PetscErrorCode    ierr;
   PetscMPIInt       rank;
   ierr = PetscInitialize(&argc,&argv,(char*)0,help);CHKERRQ(ierr);
   comm = PETSC_COMM_WORLD;
   ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);

   // Problem parameters
   AppCtx user;
   ierr = PetscStrcpy(user.outputBasename, "sol"); CHKERRQ(ierr);
   user.vtkInterval = 10;

   // Options
   char filename[PETSC_MAX_PATH_LEN] = "./tri.msh";
   ierr = PetscOptionsBegin(comm,NULL,"Options","");CHKERRQ(ierr);
   ierr = PetscOptionsString("-mesh","Gmsh filename to read","",filename,filename,sizeof(filename),NULL);CHKERRQ(ierr);
   ierr = PetscOptionsEnd();CHKERRQ(ierr);

   // Create the mesh
   DM        dm,dmDist;
   PetscInt  overlap=1;
   ierr = DMPlexCreateGmshFromFile(comm,filename,PETSC_TRUE,&dm);CHKERRQ(ierr);

   // Tell the DM how degrees of freedom interact
   ierr = DMSetBasicAdjacency(dm, PETSC_TRUE, PETSC_FALSE); CHKERRQ(ierr);

   // Distribute the mesh
   ierr = DMPlexDistribute(dm,overlap,NULL,&dmDist);CHKERRQ(ierr);
   if (dmDist) { ierr = DMDestroy(&dm);CHKERRQ(ierr); dm = dmDist; }
   ierr = DMSetFromOptions(dm);CHKERRQ(ierr);
   ierr = DMSetUp(dm); CHKERRQ(ierr);
   ierr = DMView(dm, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
   ierr = DMGetDimension(dm,&(user.dim));CHKERRQ(ierr);

   // Setup the section, 1 dof per cell
   PetscSection sec;
   PetscInt     p,pStart,pEnd;
   ierr = PetscSectionCreate(PetscObjectComm((PetscObject)dm),&sec);CHKERRQ(ierr);
   ierr = PetscSectionSetNumFields(sec,1);CHKERRQ(ierr);
   ierr = PetscSectionSetFieldName(sec,0,"u");CHKERRQ(ierr);
   ierr = PetscSectionSetFieldComponents(sec,0,1);CHKERRQ(ierr);
   ierr = DMPlexGetChart(dm,&pStart,&pEnd);CHKERRQ(ierr);
   ierr = PetscSectionSetChart(sec,pStart,pEnd);CHKERRQ(ierr);
   ierr = DMPlexGetHeightStratum(dm,0,&pStart,&pEnd);CHKERRQ(ierr);
   for(p=pStart;p<pEnd;p++)
   {
      ierr = PetscSectionSetFieldDof(sec,p,0,1); CHKERRQ(ierr);
      ierr = PetscSectionSetDof(sec,p,1); CHKERRQ(ierr);
   }
   ierr = PetscSectionSetUp(sec);CHKERRQ(ierr);
   ierr = DMSetSection(dm,sec);CHKERRQ(ierr);
   ierr = PetscSectionDestroy(&sec);CHKERRQ(ierr);

   // Setup problem info
   ierr = MeshInfoCreate(dm,&(user.info));CHKERRQ(ierr);

   // Create a vec for the initial condition
   Vec U;
   ierr = DMCreateGlobalVector(dm,&U);CHKERRQ(ierr);
   // Set initial condition
   ierr = SetIC(dm,&user,U); CHKERRQ(ierr);
   ierr = PetscObjectSetName((PetscObject)U,"u");CHKERRQ(ierr);

   // Create time stepping and solve
   TS  ts;
   ierr = TSCreate(comm,&ts);CHKERRQ(ierr);
   ierr = TSSetProblemType(ts,TS_NONLINEAR); CHKERRQ(ierr);
   ierr = TSSetType(ts,TSEULER);CHKERRQ(ierr);
   ierr = TSSetRHSFunction(ts,NULL,ComputeResidual,&user);CHKERRQ(ierr);
   ierr = TSMonitorSet(ts,Monitor,&user,NULL); CHKERRQ(ierr);
   ierr = TSSetDM(ts,dm);CHKERRQ(ierr);
   ierr = TSSetSolution(ts,U);CHKERRQ(ierr);
   ierr = TSSetMaxSteps(ts,1000000);CHKERRQ(ierr);
   ierr = TSSetMaxTime(ts,1.0);CHKERRQ(ierr);
   ierr = TSSetTimeStep(ts,0.001); CHKERRQ(ierr);
   ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_STEPOVER);CHKERRQ(ierr);
   ierr = TSSetFromOptions(ts); CHKERRQ(ierr);
   ierr = TSSetUp(ts);CHKERRQ(ierr);
   ierr = TSSolve(ts,U);CHKERRQ(ierr);

   // Cleanup
   ierr = MeshInfoDestroy(&(user.info));CHKERRQ(ierr);
   ierr = TSDestroy(&ts);CHKERRQ(ierr);
   ierr = VecDestroy(&U);CHKERRQ(ierr);
   ierr = DMDestroy(&dm);CHKERRQ(ierr);
   ierr = PetscFinalize();CHKERRQ(ierr);
   return(0);
}
