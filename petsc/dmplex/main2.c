// Run as: ./main file.msh
static char help[] = "Testing DMPlex.\n\n";

#include "petscdm.h"
#include "petscdmplex.h"
#include "petscds.h"

static PetscErrorCode OutputVTK(DM dm, const char *filename, PetscViewer *viewer)
{
   PetscViewerCreate(PetscObjectComm((PetscObject)dm), viewer);
   PetscViewerSetType(*viewer, PETSCVIEWERVTK);
   PetscViewerFileSetName(*viewer, filename);
   return(0);
}

PetscErrorCode CreatePartitionVec(DM dm, DM *dmCell, Vec *partition)
{
   PetscSF        sfPoint;
   PetscSection   coordSection;
   Vec            coordinates;
   PetscSection   sectionCell;
   PetscScalar    *part;
   PetscInt       cStart, cEnd, c;
   PetscMPIInt    rank;

   DMGetCoordinateSection(dm, &coordSection);
   DMGetCoordinatesLocal(dm, &coordinates);
   DMClone(dm, dmCell);
   DMGetPointSF(dm, &sfPoint);
   DMSetPointSF(*dmCell, sfPoint);
   DMSetCoordinateSection(*dmCell, PETSC_DETERMINE, coordSection);
   DMSetCoordinatesLocal(*dmCell, coordinates);
   MPI_Comm_rank(PetscObjectComm((PetscObject)dm), &rank);
   PetscSectionCreate(PetscObjectComm((PetscObject)dm), &sectionCell);
   DMPlexGetHeightStratum(*dmCell, 0, &cStart, &cEnd);
   PetscSectionSetChart(sectionCell, cStart, cEnd);
   for (c = cStart; c < cEnd; ++c)
   {
       PetscSectionSetDof(sectionCell, c, 1);
   }
   PetscSectionSetUp(sectionCell);
   DMSetSection(*dmCell, sectionCell);
   PetscSectionDestroy(&sectionCell);
   DMCreateLocalVector(*dmCell, partition);
   PetscObjectSetName((PetscObject)*partition, "partition");
   VecGetArray(*partition, &part);
   for (c = cStart; c < cEnd; ++c)
   {
       PetscScalar *p;
       DMPlexPointLocalRef(*dmCell, c, part, &p);
       p[0] = rank;
   }
   VecRestoreArray(*partition, &part);
   return(0);
} 

int main(int argc, char *argv[])
{
   PetscErrorCode ierr;
   DM             dm;
   PetscBool      interpolate = PETSC_TRUE;
   PetscInt       overlap = 1;
   PetscInt       cStart, cEnd, cEndInt;
   PetscMPIInt    rank;

   ierr = PetscInitialize(&argc, &argv, (char*)0, help); CHKERRQ(ierr);
   ierr = DMPlexCreateGmshFromFile(MPI_COMM_WORLD, argv[1], interpolate, &dm); CHKERRQ(ierr);

   {
      DM dmDist;
      //ierr = DMSetBasicAdjacency(dm, PETSC_TRUE, PETSC_FALSE); CHKERRQ(ierr);
      ierr = DMPlexDistribute(dm, overlap, NULL, &dmDist); CHKERRQ(ierr);
      if (dmDist)
      {
         ierr = DMDestroy(&dm); CHKERRQ(ierr);
         dm = dmDist;
      }
   }
   ierr = DMSetFromOptions(dm); CHKERRQ(ierr);

   /*
   {
      DM gdm;
      ierr = DMPlexConstructGhostCells(dm, NULL, NULL, &gdm); CHKERRQ(ierr);
      ierr = DMDestroy(&dm); CHKERRQ(ierr);
      dm   = gdm;
      ierr = DMViewFromOptions(dm, NULL, "-dm_view"); CHKERRQ(ierr);
   }
   */

   ierr = DMPlexGetDepthStratum(dm, 2, &cStart, &cEnd); CHKERRQ(ierr);
   ierr = DMPlexGetHybridBounds(dm, &cEndInt, NULL, NULL, NULL); CHKERRQ(ierr);
   ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
   printf("rank = %d, cStart, cEnd, cEndInt, cEnd-cStart = %d %d %d %d\n",
          rank, cStart, cEnd, cEndInt, cEnd-cStart);

   {
      DM          dmCell;
      Vec         partition;
      PetscViewer viewer;
      ierr = CreatePartitionVec(dm, &dmCell, &partition); CHKERRQ(ierr);
      ierr = OutputVTK(dmCell, "partition.vtk", &viewer); CHKERRQ(ierr);
      ierr = VecView(partition, viewer); CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
      ierr = VecDestroy(&partition); CHKERRQ(ierr);
      ierr = DMDestroy(&dmCell); CHKERRQ(ierr);
   }

   ierr = PetscFinalize(); CHKERRQ(ierr);
}
