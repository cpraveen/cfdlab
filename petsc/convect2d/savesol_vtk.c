//------------------------------------------------------------------------------
PetscErrorCode savesol_vtk(int *c, DM da, Vec ug)
{
   PetscErrorCode ierr;
   char           filename[32] = "sol";
   PetscMPIInt    rank;
   PetscInt       i, j, ibeg, jbeg, nlocx, nlocy;
   FILE           *fp;
   Vec            ul;
   PetscScalar    **u;

   ierr = DMGetLocalVector(da, &ul); CHKERRQ(ierr);
   ierr = DMGlobalToLocalBegin(da, ug, INSERT_VALUES, ul); CHKERRQ(ierr);
   ierr = DMGlobalToLocalEnd(da, ug, INSERT_VALUES, ul); CHKERRQ(ierr);
   ierr = DMDAVecGetArray(da, ul, &u); CHKERRQ(ierr);
   ierr = DMDAGetCorners(da, &ibeg, &jbeg, 0, &nlocx, &nlocy, 0); CHKERRQ(ierr);

   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
   sprintf(filename, "sol-%03d-%03d.vtr", rank, *c);
   fp = fopen(filename,"w");
   fprintf(fp, "# vtk DataFile Version 3.0\n");
   fprintf(fp, "Sample rectilinear grid\n");
   fprintf(fp, "ASCII\n");
   fprintf(fp, "DATASET RECTILINEAR_GRID\n");
   fprintf(fp, "DIMENSIONS %d %d %d\n", nlocx+1, nlocy+1, 1);
   fprintf(fp, "X_COORDINATES %d float\n", nlocx+1);
   for(i=ibeg; i<ibeg+nlocx+1; ++i)
   {
      PetscReal x = xmin + i*dx + 0.5*dx;
      fprintf(fp, "%e ", x);
   }
   fprintf(fp, "\n");
   fprintf(fp, "Y_COORDINATES %d float\n", nlocy+1);
   for(j=jbeg; j<jbeg+nlocy+1; ++j)
   {
      PetscReal y = ymin + j*dy + 0.5*dy;
      fprintf(fp, "%e ", y);
   }
   fprintf(fp, "\n");
   fprintf(fp, "Z_COORDINATES 1 float\n");
   fprintf(fp, "0.0\n");

   fprintf(fp, "POINT_DATA %d\n", (nlocx+1)*(nlocy+1)); 
   fprintf(fp, "SCALARS scalars float\n");
   fprintf(fp, "LOOKUP_TABLE default\n");

   for(j=jbeg; j<jbeg+nlocy+1; ++j)
      for(i=ibeg; i<ibeg+nlocx+1; ++i)
      {
         fprintf(fp,"%e ",u[j][i]);
      }
   fclose(fp);

   ierr = DMDAVecRestoreArray(da, ul, &u); CHKERRQ(ierr);
   ierr = DMRestoreLocalVector(da, &ul); CHKERRQ(ierr);

   ++(*c);
   return(0);
}
