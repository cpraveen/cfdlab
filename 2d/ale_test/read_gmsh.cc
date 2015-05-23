#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include <cstdlib>
#include "grid.h"

using namespace std;

//------------------------------------------------------------------------------
// Read a grid in gmsh format
//------------------------------------------------------------------------------
void Grid::read_gmsh (string grid_file)
{
   unsigned int i, j, count, elem_type, ntags, tag, n_elem;
   string line;

   cout << "Reading gmsh grid file " << grid_file << endl;

   ifstream file;
   file.open (grid_file.c_str());
   assert ( file.is_open() );

   // Skip some lines
   file >> line;
   file >> line;
   file >> line;
   file >> line;
   file >> line;
   file >> line;

   // Read vertices
   file >> n_vertex;
   assert (n_vertex > 0);
   vertex.resize (n_vertex);

   double z;
   for(i=0; i<n_vertex; ++i)
      file >> count
           >> vertex[i].coord.x
           >> vertex[i].coord.y
           >> z;

   file >> line;
   file >> line;

   // Read cells
   file >> n_elem;
   assert (n_elem > 0);

   for(i=0; i<n_elem; ++i)
   {
      file >> count
           >> elem_type
           >> ntags;

      // Some gmsh files have 2 and some have 3 tags
      assert( ntags==2 || ntags==3 );

      if(elem_type == 1) // Line face
      {
         face.resize (n_face+1);
         file >> face[n_face].type; // First tag is face type
         // Dummy tags
         for(unsigned int p=1; p<ntags; ++p)
            file >> tag;
         file >> face[n_face].vertex[0] 
              >> face[n_face].vertex[1];
         ++n_face;
      }
      else if(elem_type == 2) // Triangular cell
      {
         cell.resize (n_cell+1);
         for(unsigned int p = 0; p<ntags; ++p)
               file >> tag ; // Dummy tags

         file >> cell[n_cell].vertex[0] 
              >> cell[n_cell].vertex[1] 
              >> cell[n_cell].vertex[2];
         ++n_cell;
      }
      else
      {
         cout << "Unknown element type !!!" << endl;
         cout << "   Element type =" << elem_type << endl;
         exit (0);
      }
   }

   file.close ();

   cout << "No of vertices = " << n_vertex << endl;
   cout << "No of cells    = " << n_cell << endl;

   // Node number must start from 0
   // but gmsh starts from 1. So we reduce count by one
   // Vertices forming the cell
   for(i=0; i<n_cell; ++i)
      for(j=0; j<3; ++j)
         cell[i].vertex[j]--;

   // Vertices forming the face
   for(i=0; i<n_face; ++i)
      for(j=0; j<2; ++j)
         face[i].vertex[j]--;

   for(i=0; i<n_cell; ++i)
   {
      for(j=0; j<3; ++j)
      {
         int n = cell[i].vertex[j];
         cell[i].centroid += vertex[n].coord;
      }
      cell[i].centroid *= (1.0/3.0);
   }

}
