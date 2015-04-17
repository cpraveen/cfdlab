/* Convert a Cobalt mesh file to su2 format. Cobalt mesh must contain only
 * triangles and tetrahedra. 
 * Written by Praveen. C, http://math.tifrbng.res.in/~praveen 
 * E.g.:
 *     mesh    = mesh file
 *     mesh.bc = bc information file
 *     To run, do
 *         ./cobalt_to_su2 mesh
 *     It should create a file mesh.su2
*/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#define NFACE_TYPE_MAX  100

int fourth_vert(int *cvert, int *fvert);

typedef struct _FACE_ {
   int vert[3];
   int lcell;
   int rcell;
} FACE;

typedef struct _CELL_ {
   int vert[4];
   int filled;
} CELL;

int main(int argc, char *argv[]){
   FILE *fp;
   FILE *fo;
   int i, j;
   int lcell, rcell;
   int idummy;
   int npoints,nfaces,ncells,nmark;
   int nbfaces;
   int face_id;
   int facecount[NFACE_TYPE_MAX];
   int *nmark_face[NFACE_TYPE_MAX];
   int count;
   double *x, *y, *z;
   FACE *face;
   CELL *cell;
   char bc_names[NFACE_TYPE_MAX][1024];
   char bcfile[1024];
   FILE *fbc;
   char dummy[1024];
   char su2_file[1024];

   // Read bc names from *.bc file
   strcpy(bcfile, argv[1]);
   strcat(bcfile,".bc");
   printf("Reading bc info from %s\n", bcfile);
   fbc = fopen(bcfile, "r");
   // Ignore first four lines
   fgets(dummy, 1024, fbc);
   fgets(dummy, 1024, fbc);
   fgets(dummy, 1024, fbc);
   fgets(dummy, 1024, fbc);
   while(feof(fbc)==0)
   {
      fscanf(fbc, "%d", &i);
      fscanf(fbc, "%s\n", bc_names[i]);
      //printf("%d   %s\n", i, bc_names[i]);
      fgets(dummy, 1024, fbc);
      fgets(dummy, 1024, fbc);
      fgets(dummy, 1024, fbc);
      fgets(dummy, 1024, fbc);
   }
   fclose(fbc);

   for(i=0; i<NFACE_TYPE_MAX; i++) facecount[i]=0;

   printf("Reading grid file %s\n", argv[1]);
   fp = fopen(argv[1],"r");
   if(fp==NULL){
      printf("Could not open cobalt file\n");
      exit(0);
   }
   fscanf(fp,"%d%d%d",&idummy,&idummy,&idummy);
   fscanf(fp,"%d%d%d%d%d",&npoints,&nfaces,&ncells,&idummy,&idummy);
   printf("npoints = %d\n",npoints);
   printf("nfaces  = %d\n",nfaces);
   printf("ncells  = %d\n",ncells);

   x = (double*)malloc(npoints*sizeof(double));
   y = (double*)malloc(npoints*sizeof(double));
   z = (double*)malloc(npoints*sizeof(double));

   for(i=0; i<npoints; i++){
      fscanf(fp,"%lf%lf%lf",&x[i],&y[i],&z[i]);
   }

   face = (FACE*)malloc(nfaces*sizeof(FACE));
   for(i=0; i<nfaces; i++)
   {
      fscanf(fp,"%d%d%d%d%d%d",&idummy,&face[i].vert[0],&face[i].vert[1],
            &face[i].vert[2],&face[i].lcell,&face[i].rcell);
      if(idummy != 3)
      {
         printf("Found a non-triangular face with %d vertices\n", idummy);
         exit(0);
      }
   }
   fclose(fp);

   // count number of boundary faces
   nbfaces = 0;
   for(i=0; i<nfaces; i++)
      if(face[i].rcell < 0){
         assert(-face[i].rcell < NFACE_TYPE_MAX); // increase NFACE_TYPE_MAX
         ++facecount[-face[i].rcell];
         ++nbfaces;
      }
   printf("No. of boundary faces = %d\n",nbfaces);
   nmark = 0; // number of boundary face types
   
   for(i=0; i<NFACE_TYPE_MAX; i++)
      if(facecount[i] > 0)
      {
         printf("%5d %10d %s\n", i, facecount[i], bc_names[i]);
         ++nmark;
         nmark_face[i] = (int*)malloc(facecount[i]*sizeof(int));
      }
   printf("Number of face markers = %d\n", nmark);

   // reset to zero
   for(i=0; i<NFACE_TYPE_MAX; i++) facecount[i]=0;

   for(i=0; i<nfaces; i++)
   {
      if(face[i].rcell < 0)
      {
         nmark_face[-face[i].rcell][facecount[-face[i].rcell]] = i;
         ++facecount[-face[i].rcell];
      }
   }


   // create cell information
   cell = (CELL*)malloc(ncells*sizeof(CELL));
   for(i=0; i<ncells; i++) cell[i].filled = 0;

   for(i=0; i<nfaces; i++){

      lcell = face[i].lcell - 1;
      if(cell[lcell].filled==0){
         cell[lcell].vert[0] = face[i].vert[0];
         cell[lcell].vert[1] = face[i].vert[1];
         cell[lcell].vert[2] = face[i].vert[2];
         cell[lcell].filled = 1;
      }
      else if(cell[lcell].filled==1){
         cell[lcell].vert[3] = fourth_vertex(cell[lcell].vert,face[i].vert);
         cell[lcell].filled = 2;
      }

      rcell = face[i].rcell - 1;
      if(rcell+1 > 0 && cell[rcell].filled==0){
         cell[rcell].vert[0] = face[i].vert[2]; // note opposite order
         cell[rcell].vert[1] = face[i].vert[1];
         cell[rcell].vert[2] = face[i].vert[0];
         cell[rcell].filled = 1;
      }
      else if(rcell+1 > 0 && cell[rcell].filled==1){
         cell[rcell].vert[3] = fourth_vertex(cell[rcell].vert,face[i].vert);
         cell[rcell].filled = 2;
      }

   }


   strcpy(su2_file, argv[1]);
   strcat(su2_file, ".su2");
   printf("Writing su2 mesh file into %s\n", su2_file);
   fo = fopen(su2_file, "w");
   fprintf(fo,"NDIME= 3\n");

   // Save tetrahedra
   fprintf(fo,"NELEM= %d\n", ncells);
   for(i=0; i<ncells; i++)
   {
      fprintf(fo, "10  "); // VTK type for tetrahedra
      for(j=0; j<4; j++)
         fprintf(fo, "%d ", cell[i].vert[j]-1);
      fprintf(fo, "%d\n", i);
   }

   // Save point coordinates
   fprintf(fo,"NPOIN= %d\n", npoints);
   for(i=0; i<npoints; i++)
      fprintf(fo, "%lf %lf %lf %d\n", x[i], y[i], z[i], i);

   fprintf(fo,"NMARK= %d\n", nmark);
   // Save boundary faces
   for(i=0; i<NFACE_TYPE_MAX; i++)
      if(facecount[i] > 0)
      {
         //fprintf(fo,"MARKER_TAG= %d\n", i);
         fprintf(fo,"MARKER_TAG= %s\n", bc_names[i]);
         fprintf(fo,"MARKER_ELEMS= %d\n", facecount[i]);
         for(j=0; j<facecount[i]; ++j)
         {
            face_id = nmark_face[i][j];
            fprintf(fo, "5  "); // VTK type for triangle
            fprintf(fo, "%d  %d  %d\n", face[face_id].vert[0]-1,
                                        face[face_id].vert[1]-1,
                                        face[face_id].vert[2]-1);
         }
      }
   fprintf(fo, "NCHUNK= 0\n");

   fclose(fo);
   free(face);
   free(cell);

}

int fourth_vertex(int *cvert, int *fvert){
   int i;

   for(i=0; i<3; i++)
      if(fvert[i] != cvert[0] && fvert[i] != cvert[1] && fvert[i] != cvert[2])
         return fvert[i];

   printf("fourth_vertex: fatal error, fourth vertex not found");
}
