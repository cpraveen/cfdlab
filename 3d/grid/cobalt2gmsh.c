/* Convert a Cobalt mesh file to gmsh format. Cobalt mesh must contain only
 * triangles and tetrahedra. 
 * Written by Praveen. C, http://math.tifrbng.res.in/~praveen 
*/
#include <stdio.h>
#include <stdlib.h>

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
   int i;
   int lcell, rcell;
   int idummy;
   int npoints,nfaces,ncells;
   int nbfaces;
   int facecount[100];
   int count;
   double x, y, z;
   FACE *face;
   CELL *cell;
   double wing_span = 29.4; // meter
   double wing_area = 61.0; // meter^2
   double mean_chord = wing_area/wing_span;

   printf("wing span  = %e meters\n",  wing_span);
   printf("wing area  = %e meter^2\n", wing_area);
   printf("mean chord = %e meters\n",  mean_chord);

   for(i=0; i<100; i++) facecount[i]=0;

   fo = fopen("out.msh","w");
   fprintf(fo,"$MeshFormat\n");
   fprintf(fo,"2.0   0    8\n");
   fprintf(fo,"$EndMeshFormat\n");

   fp = fopen(argv[1],"r");
   if(fp==NULL){
      printf("Could not open cobalt file\n");
      exit(0);
   }
   fscanf(fp,"%d%d%d",&idummy,&idummy,&idummy);
   fscanf(fp,"%d%d%d%d%d",&npoints,&nfaces,&ncells,&idummy,&idummy);
   printf("%d %d %d\n",npoints,nfaces,ncells);

   fprintf(fo,"$Nodes\n");
   fprintf(fo,"%d\n",npoints);

   printf("Coordinates will be scaled ...\n");

   // convert from mm to meters
   // make mean chord = 1
   for(i=0; i<npoints; i++){
      fscanf(fp,"%lf%lf%lf",&x,&y,&z);
      x = x/(1000.0*mean_chord);
      y = y/(1000.0*mean_chord);
      z = z/(1000.0*mean_chord);
      fprintf(fo,"%d %lf %lf %lf\n",i+1,x,y,z);
   }
   fprintf(fo,"$EndNodes\n");

   face = (FACE*)malloc(nfaces*sizeof(FACE));
   for(i=0; i<nfaces; i++){
      fscanf(fp,"%d%d%d%d%d%d",&idummy,&face[i].vert[0],&face[i].vert[1],
            &face[i].vert[2],&face[i].lcell,&face[i].rcell);
   }
   fclose(fp);

   // count number of boundary faces
   nbfaces = 0;
   for(i=0; i<nfaces; i++)
      if(face[i].rcell < 0){
         ++facecount[-face[i].rcell];
         ++nbfaces;
      }
   printf("No. of boundary faces = %d\n",nbfaces);
   for(i=0; i<100; i++)
      if(facecount[i] > 0)
         printf("%5d %10d\n", i, facecount[i]);

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


   // write out cells
   fprintf(fo,"$Elements\n");
   fprintf(fo,"%d\n",ncells+nbfaces);
   count=0;
   for(i=0; i<nfaces; i++)
      if(face[i].rcell < 0)
         fprintf(fo,"%d %d %d %d %d %d %d %d %d\n",++count,2,3,
               -face[i].rcell,0,0,face[i].vert[0],face[i].vert[1],
               face[i].vert[2]);
   for(i=0; i<ncells; i++)
      fprintf(fo,"%d %d %d %d %d %d %d %d %d %d\n",++count,4,3,
            0,0,0,cell[i].vert[0],cell[i].vert[1],
            cell[i].vert[2],cell[i].vert[3]);
   fprintf(fo,"$EndElements\n");

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
