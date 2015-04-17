/* Convert cobalt mesh file to su2 format
 * Author: Praveen. C
 *         http://math.tifrbng.res.in/~praveen
 *
 * Elements supported: tetrahedron, wedge, pyramid, triangle, quadrilateral
 * Usage: Compile the code
 *        c++ main.cc -o cobalt2su
 * You need the two cobalt files, e.g
 *        foo
 *        foo.bc
 * Run the converter as:
 *        cobalt2su foo
 * This creates the file foo.su2 and a vtk file mesh.vtk
*/
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <fstream>
#include <cstring>

#define LEFT    0
#define RIGHT   1

#define TRI  0
#define QUAD 1

// VTK element types
#define TETRA   10 // tetrahedron
#define WEDGE   13 // prism with triangular base, 
                   // has 2 triangles, 3 quadrilateral
#define PYRAMID 14 // pyramid: 4 triangle, 1 quadrilateral

#define NFACE_TYPE_MAX   1000

using namespace std;

struct Face
{
   int type;
   vector<int> vertex;
   int lcell, rcell;
};

struct Cell
{
   int type;
   vector<int> face_id;
   vector<int> face_pos;
   vector<int> vertex;
};

void find_tetra_vertex(Cell& cell, const vector<Face>& face);
void find_wedge_vertex(Cell& cell, const vector<Face>& face);
void find_pyramid_vertex(Cell& cell, const vector<Face>& face);
int find_in_v1_notin_v2(const vector<int>& v1, const vector<int>& v2);

int main(int argc, char *argv[])
{
   unsigned int ncell, nface, nvertex;
   int dummy;

   ifstream fc;
   fc.open(argv[1]);

   fc >> dummy >> dummy >> dummy;
   fc >> nvertex >> nface >> ncell;
   cout << "ncell = " << ncell 
        << " nface = " << nface 
        << " nvertex = " << nvertex
        << "\n";
   fc >> dummy >> dummy;

   cout << "Reading coordinates ...";
   vector<double> x(nvertex), y(nvertex), z(nvertex);
   double xmin, ymin, zmin, xmax, ymax, zmax;
   xmin = ymin = zmin = +1e20;
   xmax = ymax = zmax = -1e20;
   for(unsigned int i=0; i<nvertex; ++i)
   {
      fc >> x[i] >> y[i] >> z[i];
      xmin = min(xmin, x[i]);
      xmax = max(xmax, x[i]);
      ymin = min(ymin, y[i]);
      ymax = max(ymax, y[i]);
      zmin = min(zmin, z[i]);
      zmax = max(zmax, z[i]);
   }
   cout << "Done\n";
   cout << "x min, max = " << xmin << "  " << xmax << endl;
   cout << "y min, max = " << ymin << "  " << ymax << endl;
   cout << "z min, max = " << zmin << "  " << zmax << endl;

   cout << "Reading faces ...";
   vector<Face> face(nface);
   int n_tri_face = 0;
   int n_quad_face = 0;
   int min_vertex_no=100000000;
   int min_cell_no=100000000;
   int max_vertex_no=-100000000;
   int max_cell_no=-100000000;
   for(unsigned int i=0; i<nface; ++i)
   {
      unsigned int n_face_points;
      fc >> n_face_points;
      assert(n_face_points==3 || n_face_points==4);
      face[i].vertex.resize(n_face_points);
      //for(unsigned int j=0; j<n_face_points; ++j)
      for(int j=n_face_points-1; j>=0; --j)
      {
         int p;
         fc >> p;
         min_vertex_no = min(min_vertex_no, p);
         max_vertex_no = max(max_vertex_no, p);
         --p;
         //face[i].vertex.push_back(p);
         face[i].vertex[j] = p;
      }
      fc >> face[i].lcell >> face[i].rcell;
      assert(face[i].lcell > 0);
      min_cell_no = min(min_cell_no, face[i].lcell);
      max_cell_no = max(max_cell_no, face[i].lcell);
      --face[i].lcell;
      assert(face[i].rcell != 0);
      if(face[i].rcell > 0)
      {
         min_cell_no = min(min_cell_no, face[i].rcell);
         max_cell_no = max(max_cell_no, face[i].rcell);
         --face[i].rcell;
      }
      if(n_face_points == 3)
      {
         face[i].type = TRI;
         ++n_tri_face;
      }
      else if(n_face_points == 4)
      {
         face[i].type = QUAD;
         ++n_quad_face;
      }
      else
      {
         cout << "Face has " << n_face_points << " points\n";
         exit(0);
      }
   }
   cout << "Done\n";
   fc.close();
   cout << "Minimum vertex number = " << min_vertex_no << "\n";
   cout << "Maximum vertex number = " << max_vertex_no << "\n";
   cout << "Minimum cell   number = " << min_cell_no << "\n";
   cout << "Maximum cell   number = " << max_cell_no << "\n";

   cout << "Number of triangular    faces = " << n_tri_face << "\n";
   cout << "Number of quadrilateral faces = " << n_quad_face << "\n";

   cout << "Adding faces to cells ...";
   vector<Cell> cell(ncell);
   for(unsigned int i=0; i<nface; ++i)
   {
      cell[face[i].lcell].face_id.push_back(i);
      cell[face[i].lcell].face_pos.push_back(LEFT);
      if(face[i].rcell >= 0)
      {
         cell[face[i].rcell].face_id.push_back(i);
         cell[face[i].rcell].face_pos.push_back(RIGHT);
      }
   }
   cout << "Done\n";

   cout << "Finding cell type ...";
   int n_wedge=0, n_tetra=0, n_pyramid=0;
   for(unsigned int i=0; i<ncell; ++i)
   {
      int ntri=0, nquad=0;
      for(unsigned int j=0; j<cell[i].face_id.size(); ++j)
      {
         int f = cell[i].face_id[j];
         if(face[f].type == TRI)
            ++ntri;
         else if(face[f].type == QUAD)
            ++nquad;
         else
         {
            cout << "Unknown face type\n";
            exit(0);
         }
      }

      if(ntri==4 && nquad==0)
      {
         cell[i].type = TETRA;
         ++n_tetra;
         find_tetra_vertex(cell[i], face);
      }
      else if(ntri==2 && nquad==3)
      {
         cell[i].type = WEDGE;
         ++n_wedge;
         find_wedge_vertex(cell[i], face);
      }
      else if(ntri==4 && nquad==1)
      {
         cell[i].type = PYRAMID;
         ++n_pyramid;
         find_pyramid_vertex(cell[i], face);
      }
      else
      {
         cout << "Unknown cell type, cell = " << i << "\n";
         cout << "ntri  = " << ntri  << "\n";
         cout << "nquad = " << nquad << "\n";
         exit(0);
      }
   }
   cout << "Done\n";

   cout << "Number of tetrahedra = " << n_tetra   << "\n";
   cout << "Number of wedge      = " << n_wedge   << "\n";
   cout << "Number of pyramid    = " << n_pyramid << "\n";

   int minp=+100000;
   int maxp=-1;
   for(unsigned int i=0; i<ncell; ++i)
      for(unsigned int j=0; j<cell[i].vertex.size(); ++j)
      {
         minp = min(minp, cell[i].vertex[j]);
         maxp = max(maxp, cell[i].vertex[j]);
      }

    cout << "Min vertex no in cell = " << minp << "\n";
    cout << "Max vertex no in cell = " << maxp << "\n";

   int minf=+100000;
   int maxf=-1;
   for(unsigned int i=0; i<nface; ++i)
      for(unsigned int j=0; j<face[i].vertex.size(); ++j)
      {
         minf = min(minf, face[i].vertex[j]);
         maxf = max(maxf, face[i].vertex[j]);
      }
   cout << "Min vertex no in face = " << minf << "\n";
   cout << "Max vertex no in face = " << maxf << "\n";

   // Read the bc file
   char bc_names[NFACE_TYPE_MAX][1024];
   char bcfile[1024];
   FILE *fbc;
   char cdummy[1024];

   // Read bc names from *.bc file
   strcpy(bcfile, argv[1]);
   strcat(bcfile,".bc");
   printf("Reading bc info from %s\n", bcfile);
   fbc = fopen(bcfile, "r");
   // Ignore first four lines
   fgets(cdummy, 1024, fbc);
   fgets(cdummy, 1024, fbc);
   fgets(cdummy, 1024, fbc);
   fgets(cdummy, 1024, fbc);
   while(feof(fbc)==0)
   {
      int i;
      assert(i > 0);
      fscanf(fbc, "%d", &i);
      fscanf(fbc, "%s\n", bc_names[i]);
      fgets(cdummy, 1024, fbc);
      fgets(cdummy, 1024, fbc);
      fgets(cdummy, 1024, fbc);
      fgets(cdummy, 1024, fbc);
   }
   fclose(fbc);

   int facecount[NFACE_TYPE_MAX];
   for(unsigned int i=0; i<NFACE_TYPE_MAX; i++) facecount[i]=0;

   // count number of boundary faces
   int nbfaces = 0;
   for(unsigned int i=0; i<nface; i++)
      if(face[i].rcell < 0){
         assert(-face[i].rcell < NFACE_TYPE_MAX); // increase NFACE_TYPE_MAX
         ++facecount[-face[i].rcell];
         ++nbfaces;
      }
   printf("No. of boundary faces = %d\n",nbfaces);
   int nmark = 0; // number of boundary face types

   int *nmark_face[NFACE_TYPE_MAX];
   for(unsigned int i=0; i<NFACE_TYPE_MAX; i++)
      if(facecount[i] > 0)
      {
         printf("%5d %10d %s\n", i, facecount[i], bc_names[i]);
         ++nmark;
         nmark_face[i] = (int*)malloc(facecount[i]*sizeof(int));
      }
   printf("Number of face markers = %d\n", nmark);

   // reset to zero
   for(unsigned int i=0; i<NFACE_TYPE_MAX; i++) facecount[i]=0;

   int n_bd_tri = 0;
   int n_bd_quad= 0;
   for(unsigned int i=0; i<nface; i++)
   {
      if(face[i].rcell < 0)
      {
         nmark_face[-face[i].rcell][facecount[-face[i].rcell]] = i;
         ++facecount[-face[i].rcell];
         if(face[i].type == TRI)  ++n_bd_tri;
         if(face[i].type == QUAD) ++n_bd_quad;
      }
   }
   cout << "No of boundary triangles      = " << n_bd_tri  << endl;
   cout << "No of boundary quadrilaterals = " << n_bd_quad << endl;

   //----------------------------------------------------------------
   // Now save to su2 file
   //----------------------------------------------------------------
   char su2_file[1024];
   strcpy(su2_file, argv[1]);
   strcat(su2_file, ".su2");
   printf("Writing su2 mesh file into %s\n", su2_file);
   FILE *fo;
   fo = fopen(su2_file, "w");
   fprintf(fo,"NDIME= 3\n");

   // Save tetrahedra
   fprintf(fo,"NELEM= %d\n", ncell);
   for(unsigned int i=0; i<ncell; i++)
   {
      fprintf(fo, "%d  ", cell[i].type); // VTK type for tetrahedra
      for(unsigned int j=0; j<cell[i].vertex.size(); j++)
         fprintf(fo, "%d ", cell[i].vertex[j]);
      fprintf(fo, "%d\n", i);
   }

   double sfactor;
   cout << "Scaling for coordinates (x,y,z multiplied by this factor) = ";
   cin >> sfactor;
   assert(sfactor>0.0);

   // Save point coordinates
   fprintf(fo,"NPOIN= %d\n", nvertex);
   for(unsigned int i=0; i<nvertex; i++)
   {
      x[i] *= sfactor;
      y[i] *= sfactor;
      z[i] *= sfactor;
      fprintf(fo, "%24.14e %24.14e %24.14e %d\n", x[i], y[i], z[i], i);
   }

   fprintf(fo,"NMARK= %d\n", nmark);
   // Save boundary faces
   for(unsigned int i=0; i<NFACE_TYPE_MAX; i++)
      if(facecount[i] > 0)
      {
         fprintf(fo,"MARKER_TAG= %s\n", bc_names[i]);
         fprintf(fo,"MARKER_ELEMS= %d\n", facecount[i]);
         for(unsigned int j=0; j<facecount[i]; ++j)
         {
            int face_id = nmark_face[i][j];
            if(face[face_id].type == TRI)
            {
               fprintf(fo, "5  "); // VTK type for triangle
               fprintf(fo, "%d  %d  %d\n", face[face_id].vertex[0],
                                           face[face_id].vertex[1],
                                           face[face_id].vertex[2]);
            }
            else if(face[face_id].type == QUAD)
            {
               fprintf(fo, "9  "); // VTK type for triangle
               fprintf(fo, "%d  %d  %d  %d\n", face[face_id].vertex[0],
                                               face[face_id].vertex[1],
                                               face[face_id].vertex[2],
                                               face[face_id].vertex[3]);
            }
            else
            {
               cout << "Unknown face type\n";
               exit(0);
            }
         }
      }
   fprintf(fo, "NCHUNK= 0\n");

   //----------------------------------------------------------------
   // write vtk file
   //----------------------------------------------------------------
   ofstream vtk;
   vtk.open("mesh.vtk");
   vtk << "# vtk DataFile Version 3.0" << endl;
   vtk << "su2 converted mesh" << endl;
   vtk << "ASCII" << endl;
   vtk << "DATASET UNSTRUCTURED_GRID" << endl;
   vtk << "POINTS  " << nvertex << "  float" << endl;

   for(unsigned int i=0; i<nvertex; i++)
      vtk << x[i] << " " << y[i] << " " << z[i] << endl;

   int ncell_data = 5 * n_tetra + 6 * n_pyramid + 7 * n_wedge;
   vtk << "CELLS  " << ncell << "  " << ncell_data << endl;

   for(unsigned int i=0; i<ncell; i++)
   {
      vtk << cell[i].vertex.size() << " ";
      for(unsigned int j=0; j<cell[i].vertex.size(); ++j)
         vtk << cell[i].vertex[j] << " ";
      vtk << endl;
   }

   vtk << "CELL_TYPES " << ncell << endl;
   for(unsigned int i=0; i<ncell; i++)
      vtk << cell[i].type << endl;

   vtk.close();
   cout << "Wrote vtk file mesh.vtk\n";

   return 0;
}

//-----------------------------------------------------------------------------
// Find vertices belonging to tetrahedral cell
//-----------------------------------------------------------------------------
void find_tetra_vertex(Cell& cell, const vector<Face>& face)
{
   // Take first face and add its vertices to cell
   // Stand inside cell and look at vertices of face f
   int f = cell.face_id[0];
   if(cell.face_pos[0] == LEFT)
   {
      // Numbering of vertices is ccw; we add in same order
      for(unsigned int j=0; j<=2; ++j)
         cell.vertex.push_back(face[f].vertex[j]);
   }
   else
   {
      // Numbering of vertices is cw; we add in reverse order
      for(int j=2; j>=0; --j)
         cell.vertex.push_back(face[f].vertex[j]);
   }

   // Collect vertices belonging to other 3 faces into v
   vector<int> v;
   for(unsigned int j=1; j<4; ++j)
   {
      int ff = cell.face_id[j];
      for(unsigned int k=0; k<3; ++k)
         v.push_back(face[ff].vertex[k]);
   }

   // Find fourth vertex
   cell.vertex.push_back( find_in_v1_notin_v2(v, cell.vertex) );

}

//-----------------------------------------------------------------------------
// Find vertices belonging to wedge cell
//-----------------------------------------------------------------------------
void find_wedge_vertex(Cell& cell, const vector<Face>& face)
{
   vector<int> tri_face;
   vector<int> quad_face;
   for(unsigned int i=0; i<5; ++i)
   {
      if(face[cell.face_id[i]].type == TRI) 
         tri_face.push_back( i );
      else if(face[cell.face_id[i]].type == QUAD) 
         quad_face.push_back( i );
   }
   assert( tri_face.size()  == 2);
   assert( quad_face.size() == 3);

   // Take first triangular face
   const int f    = cell.face_id[ tri_face[0] ];
   const int floc = tri_face[0];
   if(cell.face_pos[floc] == LEFT)
   {
      // Numbering of vertices is cw; we add in reverse order
      for(int j=2; j>=0; --j)
         cell.vertex.push_back(face[f].vertex[j]);
   }
   else
   {
      // Numbering of vertices is ccw; we add in same order
      for(unsigned int j=0; j<=2; ++j)
         cell.vertex.push_back(face[f].vertex[j]);
   }

   // Take second triangular face
   const int f2   = cell.face_id[ tri_face[1] ];

   // Loop over first three points in cell
   // and find corresponding point in other triangle
   for(unsigned int k=0; k<3; ++k)
   {
      // Find point in f2 connected to p0
      const int p0 = cell.vertex[k];
      bool found = false;
      int c  = 0;
      while(!found && c<3)
      {
         int both_belong = 0;
         // loop over quad faces
         for(unsigned int i=0; i<3; ++i)
         {
            int c_p0 = 0;
            int c_c  = 0;
            int f_id = cell.face_id[ quad_face[i] ];
            for(unsigned int j=0; j<4; ++j) // loop over vertices of quad face
            {
               if(face[f_id].vertex[j] == p0) ++c_p0;
               if(face[f_id].vertex[j] == face[f2].vertex[c]) ++c_c;
            }
            if(c_p0 == 1 && c_c == 1) ++both_belong;
         }
         if(both_belong == 2)
         {
            found = true;
            ++c;
         }
         else
            ++c;
      }
      --c;
      assert(found && c<3);
      cell.vertex.push_back(face[f2].vertex[c]);
   }

}

//-----------------------------------------------------------------------------
// Find vertices belonging to pyramid cell
//-----------------------------------------------------------------------------
void find_pyramid_vertex(Cell& cell, const vector<Face>& face)
{
   vector<int> tri_face;
   vector<int> quad_face;
   for(unsigned int i=0; i<5; ++i)
   {
      if(face[cell.face_id[i]].type == TRI) 
         tri_face.push_back( i );
      else if(face[cell.face_id[i]].type == QUAD) 
         quad_face.push_back( i );
   }
   assert( tri_face.size()  == 4);
   assert( quad_face.size() == 1);

   // Take quad face
   const int f    = cell.face_id[ quad_face[0] ];
   const int floc = quad_face[0];
   if(cell.face_pos[floc] == LEFT)
   {
      // Numbering of vertices is cw; we add in reverse order
      for(int j=3; j>=0; --j)
         cell.vertex.push_back(face[f].vertex[j]);
   }
   else
   {
      // Numbering of vertices is ccw; we add in same order
      for(unsigned int j=0; j<=3; ++j)
         cell.vertex.push_back(face[f].vertex[j]);
   }

   // We have four vertices, need to find fifth
   // Take first triangular face
   const int f2   = cell.face_id[ tri_face[0] ];

   // Collect vertices belonging to f2
   vector<int> v;
   for(unsigned int k=0; k<3; ++k)
      v.push_back(face[f2].vertex[k]);

   // Find fourth vertex
   cell.vertex.push_back( find_in_v1_notin_v2(v, cell.vertex) );

}

//-----------------------------------------------------------------------------
// find element in v1 which is not in v2
// We assume that v1 contains only one element which is not in v2
//-----------------------------------------------------------------------------
int find_in_v1_notin_v2(const vector<int>& v1, const vector<int>& v2)
{
   for(unsigned int i=0; i<v1.size(); ++i)
   {
      bool status = true;
      for(unsigned int j=0; j<v2.size(); ++j)
         if(v1[i] == v2[j]) status = false;
      if(status) return v1[i];
   }
   cout << "Fatal: could not find in v1 not in v2\n";
   exit(0);
}
