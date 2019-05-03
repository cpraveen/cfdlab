// Run as: ./main file.msh
static char help[] = "Testing DMPlex.\n\n";

#include "petscdmplex.h"

int main(int argc, char *argv[])
{
   PetscErrorCode ierr;
   DM             dm;
   PetscBool      interpolate = PETSC_TRUE;
   PetscInt       dim;
   PetscInt       pStart, pEnd;
   PetscInt       vStart, vEnd;
   PetscInt       eStart, eEnd;
   PetscInt       cStart, cEnd;
   PetscSection   s;

   ierr = PetscInitialize(&argc, &argv, (char*)0, help); CHKERRQ(ierr);
   ierr = DMPlexCreateGmshFromFile(MPI_COMM_WORLD, argv[1], interpolate, &dm); CHKERRQ(ierr);

   ierr = DMGetDimension(dm, &dim); CHKERRQ(ierr);
   printf("dim = %d\n", dim);

   ierr = DMPlexGetChart(dm, &pStart, &pEnd); CHKERRQ(ierr);
   printf("chart: pstart, pend = %d %d\n", pStart, pEnd);

   // vertices
   ierr = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd); CHKERRQ(ierr);
   printf("Depth 0: vstart, vend, len = %d %d %d\n", vStart, vEnd, vEnd-vStart);
   // edges
   ierr = DMPlexGetDepthStratum(dm, 1, &eStart, &eEnd); CHKERRQ(ierr);
   printf("Depth 1: estart, eend, len = %d %d %d\n", eStart, eEnd, eEnd-eStart);
   // cells
   ierr = DMPlexGetDepthStratum(dm, 2, &cStart, &cEnd); CHKERRQ(ierr);
   printf("Depth 2: cstart, cend, len = %d %d %d\n", cStart, cEnd, cEnd-cStart);

   // cells adjacent to face
   {
      FILE * fid = fopen("face_nbr.txt","w");
      for(PetscInt e=eStart; e<eEnd; ++e)
      {
         PetscInt nbr;
         ierr = DMPlexGetSupportSize(dm, e, &nbr); CHKERRQ(ierr);

         const PetscInt *nbcells;
         ierr = DMPlexGetSupport(dm, e, &nbcells); CHKERRQ(ierr);

         if(nbr == 1) // boundary face
            fprintf(fid, "%d %d\n",e-eStart+1,nbcells[0]-cStart+1);
         else if(nbr == 2) // interior face
            fprintf(fid, "%d %d %d\n",e-eStart+1,nbcells[0]-cStart+1,nbcells[1]-cStart+1);
         else
         {
            printf("nbr is not 1 or 2\n");
            exit(0);
         }
      }
      fclose(fid);
   }

   // point coordinates
   {
      Vec coordinates;
      ierr = DMGetCoordinatesLocal(dm, &coordinates); CHKERRQ(ierr);

      const PetscScalar *coords;
      ierr = VecGetArrayRead(coordinates, &coords); CHKERRQ(ierr);

      DM dmCoord;
      ierr = DMGetCoordinateDM(dm, &dmCoord); CHKERRQ(ierr);

      FILE * fid = fopen("vertices.txt","w");
      for(PetscInt v=vStart; v<vEnd; ++v)
      {
         PetscScalar  *vertex;
         ierr = DMPlexPointLocalRead(dmCoord, v, coords, &vertex); CHKERRQ(ierr);
         fprintf(fid, "%f %f\n", vertex[0], vertex[1]);
      }
      fclose(fid);

      ierr = VecRestoreArrayRead(coordinates, &coords); CHKERRQ(ierr);
   }

   // compute cell and face geometry
   {
      Vec cellgeom, facegeom;
      ierr = DMPlexComputeGeometryFVM(dm, &cellgeom, &facegeom); CHKERRQ(ierr);

      DM dmCell;
      ierr = VecGetDM(cellgeom, &dmCell); CHKERRQ(ierr);

      const PetscScalar *cgeom;
      ierr = VecGetArrayRead(cellgeom, &cgeom); CHKERRQ(ierr);

      FILE * fid = fopen("cells.txt","w");
      for(PetscInt c=cStart; c<cEnd; ++c)
      {
         // cell properties like volume, centroid
         PetscFVCellGeom *cg;
         ierr = DMPlexPointLocalRead(dmCell, c, cgeom, &cg); CHKERRQ(ierr);
         fprintf(fid, "%d %f %f %f\n", c-cStart+1, cg->volume, 
                 cg->centroid[0], cg->centroid[1]);
      }
      fclose(fid);
      ierr = VecRestoreArrayRead(cellgeom, &cgeom); CHKERRQ(ierr);
   }

   // create section with one variable in each cell
   ierr = PetscSectionCreate(PetscObjectComm((PetscObject)dm), &s); CHKERRQ(ierr);
   ierr = PetscSectionSetChart(s, pStart, pEnd); CHKERRQ(ierr);
   for(PetscInt c=cStart; c<cEnd; ++c)
   {
      ierr = PetscSectionSetDof(s, c, 1); CHKERRQ(ierr);
   }
   ierr = PetscSectionSetUp(s); CHKERRQ(ierr);

   // create vector to store solution
   Vec lv, gv;
   ierr = DMSetDefaultSection(dm, s); CHKERRQ(ierr);
   ierr = DMGetLocalVector(dm, &lv); CHKERRQ(ierr);
   ierr = DMGetGlobalVector(dm, &gv); CHKERRQ(ierr);
}
