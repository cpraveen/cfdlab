//
//  gmsh.cpp
//  
//
//  Created by Praveen Chandrashekar on 09/09/19.
//

#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

double some_function(double x, double y)
{
   return sin(2.0*M_PI*x) * sin(2.0*M_PI*y);
}

int main()
{
   const int dim = 2;
   double *coord;
   unsigned int *cell1, *cell2, *bface;
   int *bface_type;
   unsigned int n_vertex, n_tri, n_quad, n_cell, n_bface;

   // Temporary arrays
   unsigned int n_elem;
   vector<unsigned int> elem1, elem2, elem_type, elem_phy;

   string grid_file = "cylinder.msh";
   cout << "Reading gmsh grid file " << grid_file << endl;

   ifstream file;
   file.open (grid_file.c_str());
   assert ( file.is_open() );

   string line;
   int tag;
   unsigned int count, ntags;
   double rdummy;

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
   cout << "Number of vertices = " << n_vertex << endl;
   coord = new double[dim*n_vertex];

   for(unsigned int i=0; i<n_vertex; ++i)
      file >> count
           >> coord[i*dim]
           >> coord[i*dim+1]
           >> rdummy; // 2d, dont need z coordinate

   // Skip two lines
   file >> line;
   file >> line;

   // Read cells
   file >> n_elem;
   assert (n_elem > 0);
   cout << "Numbers of elements = " << n_elem << endl;

   // for triangle and quad only
   elem1.resize(n_elem+1);
   elem_type.resize(n_elem);
   elem_phy.resize(n_elem);
   elem1[0] = 0;

   n_tri = 0;
   n_quad = 0;
   n_bface = 0;
   for(unsigned int i=0; i<n_elem; ++i)
   {
      file >> count
           >> elem_type[i]
           >> ntags;

      // Some gmsh files have 2 and some have 3 tags
      assert( ntags==2 || ntags==3 );

      if(elem_type[i] == 1) // Line face
      {
         file >> elem_phy[i]; // First tag is physical type
         // Dummy tags
         for(unsigned int p=1; p<ntags; ++p)
            file >> tag;

         elem1[i+1] = elem1[i] + 2;
         elem2.resize (elem1[i+1]);
         file >> elem2[elem1[i]]
              >> elem2[elem1[i]+1];
         ++n_bface;
      }
      else if(elem_type[i] == 2) // Triangular cell
      {
         file >> elem_phy[i]; // First tag is physical type
         for(unsigned int p = 1; p<ntags; ++p)
            file >> tag ; // Dummy tags

         elem1[i+1] = elem1[i] + 3;
         elem2.resize (elem1[i+1]);
         file >> elem2[elem1[i]]
              >> elem2[elem1[i]+1]
              >> elem2[elem1[i]+2];
         ++n_tri;
      }
      else if(elem_type[i] == 3) // Quadrilateral cell
      {
         file >> elem_phy[i]; // First tag is physical type
         for(unsigned int p = 1; p<ntags; ++p)
            file >> tag ; // Dummy tags

         elem1[i+1] = elem1[i] + 4;
         elem2.resize (elem1[i+1]);
         file >> elem2[elem1[i]]
              >> elem2[elem1[i]+1]
              >> elem2[elem1[i]+2]
              >> elem2[elem1[i]+3];
         ++n_quad;
      }
      else
      {
         cout << "Unknown element type !!!" << endl;
         cout << "   Element type =" << elem_type[i] << endl;
         exit (0);
      }
   }
   file.close();

   n_cell = n_tri + n_quad;

   cout << "Numbers of boundary faces = " << n_bface << endl;
   cout << "Number of triangles       = " << n_tri << endl;
   cout << "Number of quadrilaterals  = " << n_quad << endl;
   cout << "Total number of cells     = " << n_cell << endl;

   assert (n_cell > 0);
   assert (n_bface > 0);
   assert (elem1[n_elem] == 2*n_bface + 3*n_tri + 4*n_quad);

   // gmsh vertex numbering starts at 1, decrease by 1
   for(unsigned int i=0; i<elem1[n_elem]; ++i)
      --elem2[i];

   cell1 = new unsigned int[n_cell+1];
   cell2 = new unsigned int[3*n_tri + 4*n_quad];
   bface = new unsigned int[2*n_bface];
   bface_type = new int[n_bface];

   cell1[0] = 0;
   unsigned int bface_count = 0;
   unsigned int cell_count = 0;

   for(unsigned int i=0; i<n_elem; ++i)
   {
      if(elem_type[i] == 1) // Line face
      {
         bface_type[bface_count] = elem_phy[i];
         bface[2*bface_count] = elem2[elem1[i]];
         bface[2*bface_count+1] = elem2[elem1[i]+1];
         ++bface_count;
      }
      else if(elem_type[i] == 2 || elem_type[i] == 3) // tri and quad
      {
         cell1[cell_count+1] = cell1[cell_count] + elem1[i+1] - elem1[i];
         for(unsigned int j=0; j<elem1[i+1]-elem1[i]; ++j)
            cell2[cell1[cell_count]+j] = elem2[elem1[i]+j];
         ++cell_count;
      }
      else
      {
         cout << "Unknown element type = " << elem_type[i] << endl;
         exit(0);
      }
   }

   assert(bface_count == n_bface);
   assert(cell_count == n_cell);

   elem1.resize(0);
   elem2.resize(0);
   elem_type.resize(0);
   elem_phy.resize(0);

   // Create some data
   double *cvalue = new double[n_cell];
   double *pvalue = new double[n_vertex];

   // values at vertices
   for(unsigned int i=0; i<n_vertex; ++i)
   {
      double x = coord[i*dim];
      double y = coord[i*dim+1];
      pvalue[i] = some_function(x, y);
   }

   // values at cell centers
   for(unsigned int i=0; i<n_cell; ++i)
   {
      double x = 0;
      double y = 0;
      for(unsigned int j=cell1[i]; j<cell1[i+1]; ++j)
      {
         unsigned int v = cell2[j];
         x += coord[v*dim];
         y += coord[v*dim+1];
      }
      x /= cell1[i+1] - cell1[i];
      y /= cell1[i+1] - cell1[i];
      cvalue[i] = some_function(x, y);
   }


   // Save vtk file
   string filename = "foo.vtk";
   ofstream vtk (filename);

   vtk << "# vtk DataFile Version 3.0" << endl;
   vtk << "Test file" << endl;
   vtk << "ASCII" << endl;
   vtk << "DATASET UNSTRUCTURED_GRID" << endl;
   vtk << "POINTS  " << n_vertex << "  float" << endl;

   for(unsigned int i=0; i<n_vertex; ++i)
      vtk << coord[i*dim] << " "
          << coord[i*dim+1] << " "
          << 0.0 << endl;

   vtk << "CELLS  " << n_cell << " " << 4 * n_tri + 5 * n_quad << endl;
   for(unsigned int i=0; i<n_cell; ++i)
   {
      vtk << cell1[i+1]-cell1[i];
      for(unsigned int j=cell1[i]; j<cell1[i+1]; ++j)
         vtk << "  " << cell2[j];
      vtk << endl;
   }

   vtk << "CELL_TYPES  " << n_cell << endl;
   for(unsigned int i=0; i<n_cell; ++i)
   {
      unsigned int nv = cell1[i+1] - cell1[i];
      if(nv == 3) vtk << 5 << endl;
      else if(nv == 4) vtk << 9 << endl;
      else
      {
         cout << "Unknown nv = " << nv << endl;
         exit(0);
      }
   }

   vtk << "POINT_DATA " << n_vertex << endl;
   vtk << "SCALARS point_density float 1" << endl;
   vtk << "LOOKUP_TABLE default" << endl;
   for(unsigned int i=0; i<n_vertex; ++i)
      vtk << pvalue[i] << endl;

   vtk << "CELL_DATA " << n_cell << endl;
   vtk << "SCALARS cell_density float 1" << endl;
   vtk << "LOOKUP_TABLE default" << endl;
   for(unsigned int i=0; i<n_cell; ++i)
      vtk << cvalue[i] << endl;

   vtk.close();
   cout << "Wrote " << filename << endl;
}
