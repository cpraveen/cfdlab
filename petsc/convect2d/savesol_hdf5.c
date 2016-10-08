//------------------------------------------------------------------------------
// Saves solution into files sol000.h5, sol001.h5, etc.
// You can open them in VisIt.
//------------------------------------------------------------------------------
PetscErrorCode savesol(int *c, Vec ug)
{
   PetscErrorCode ierr;
   char           filename[32] = "sol";
   PetscViewer    viewer;
   sprintf(filename, "sol%03d.h5", *c);
   ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);  CHKERRQ(ierr);
   ierr = VecView(ug, viewer); CHKERRQ(ierr);
   ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
   ++(*c);
   return(0);
}
