#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <utility> // std::pair, std::make_pair

#include "grid.h"

using namespace std;

// constructor
Grid::Grid()
{
   has_esup = false;
   has_psup = false;
}

// desrructor
Grid::~Grid()
{
   if(!coord) delete [] coord;
   if(!cell1) delete [] cell1;
   if(!cell2) delete [] cell2;
   if(!bface) delete [] bface;
   if(!bface_type) delete [] bface_type;
   if(!esup1) delete [] esup1;
   if(!esup2) delete [] esup2;
   if(!psup1) delete [] psup1;
   if(!psup2) delete [] psup2;
}

// Read mesh from file in gmsh format.
// Currently it can only read msh2 format.
void Grid::read_gmsh(const string grid_file)
{
   // Temporary arrays
   unsigned int n_elem;
   vector<unsigned int> elem1, elem2, elem_type, elem_phy;

   cout << "Reading gmsh grid file " << grid_file << endl;

   ifstream file;
   file.open(grid_file.c_str());
   assert(file.is_open());

   string line;
   int tag;
   unsigned int count, ntags;
   double rdummy, version;

   // Skip some lines
   file >> line;
   file >> version >> tag >> tag;
   assert(version == 2.2);
   file >> line;
   file >> line;

   // Read vertices
   file >> n_vertex;
   assert(n_vertex > 0);
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
   assert(n_elem > 0);
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
      assert(ntags==2 || ntags==3);

      if(elem_type[i] == 1) // Line face
      {
         file >> elem_phy[i]; // First tag is physical type
         // Dummy tags
         for(unsigned int p=1; p<ntags; ++p)
            file >> tag;

         elem1[i+1] = elem1[i] + 2;
         elem2.resize(elem1[i+1]);
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
         elem2.resize(elem1[i+1]);
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
         elem2.resize(elem1[i+1]);
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
         exit(0);
      }
   }
   file.close();

   n_cell = n_tri + n_quad;

   cout << "Numbers of boundary faces = " << n_bface << endl;
   cout << "Number of triangles       = " << n_tri << endl;
   cout << "Number of quadrilaterals  = " << n_quad << endl;
   cout << "Total number of cells     = " << n_cell << endl;

   assert(n_cell > 0);
   assert(n_bface > 0);
   assert(elem1[n_elem] == 2*n_bface + 3*n_tri + 4*n_quad);

   // gmsh vertex numbering starts at 1, decrease by 1
   for(unsigned int i=0; i<elem1[n_elem]; ++i)
      --elem2[i];

   cell1 = new unsigned int[3*n_tri + 4*n_quad];
   cell2 = new unsigned int[n_cell+1];
   bface = new unsigned int[2*n_bface];
   bface_type = new int[n_bface];
   ctype = new int[n_cell];

   cell2[0] = 0;
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
         cell2[cell_count+1] = cell2[cell_count] + elem1[i+1] - elem1[i];
         for(unsigned int j=0; j<elem1[i+1]-elem1[i]; ++j)
            cell1[cell2[cell_count]+j] = elem2[elem1[i]+j];
         ctype[cell_count] = elem_type[i];
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

   // Delete local arrays
   elem1.resize(0);
   elem2.resize(0);
   elem_type.resize(0);
   elem_phy.resize(0);
}

// Write grid in vtk format
void Grid::write_vtk(const string grid_file)
{
   // Save vtk file
   ofstream vtk;
   vtk.open(grid_file);

   vtk << "# vtk DataFile Version 3.0" << endl;
   vtk << "Test file" << endl;
   vtk << "ASCII" << endl;
   vtk << "DATASET UNSTRUCTURED_GRID" << endl;
   vtk << "POINTS  " << n_vertex << "  float" << endl;

   for(unsigned int i=0; i<n_vertex; ++i)
   {
      auto x = get_coord(i);
      vtk << x[0] << " " << x[1] << " " << 0.0 << endl;
   }

   vtk << "CELLS  " << n_cell << " " << 4 * n_tri + 5 * n_quad << endl;
   for(unsigned int i=0; i<n_cell; ++i)
   {
      auto cell = get_cell_vertices(i);
      vtk << cell.first;
      for(unsigned int j=0; j<cell.first; ++j)
         vtk << "  " << cell.second[j];
      vtk << endl;
   }

   vtk << "CELL_TYPES  " << n_cell << endl;
   for(unsigned int i=0; i<n_cell; ++i)
   {
      if(ctype[i] == 2) vtk << 5 << endl;
      else if(ctype[i] == 3) vtk << 9 << endl;
      else
      {
         cout << "Unknown ctype = " << ctype[i] << endl;
         exit(0);
      }
   }

   vtk.close();
   cout << "Wrote " << grid_file << endl;
}

// Find cells surrounding a point
// See Lohner: section xxx
void Grid::construct_esup()
{
   if(has_esup) return;
   cout << "Constructing elements surrounding point ... ";

   esup2 = new unsigned int[n_vertex+1];
   for(unsigned int i=0; i<n_vertex+1; ++i) // Initialize esup2
   {
      esup2[i] = 0;
   }

   for(unsigned int icell=0; icell<n_cell; ++icell) // Loop over all the cells
   {
      auto cell = get_cell_vertices(icell);
      for(unsigned int ipoint=0; ipoint<cell.first; ++ipoint)  //Loop over nodes of the cell
      {
         ++esup2[cell.second[ipoint]+1]; // Count the cells surrounding the node
      }
   }

   for(unsigned int i=0; i<n_vertex; ++i) // Add the previous element to create the esup2
   {
      esup2[i+1] += esup2[i];
   }

   esup1 = new unsigned int[esup2[n_vertex]]; // Create and initialize esup1
   for(unsigned int i=0; i<esup2[n_vertex]; ++i)
   {
      esup1[i] = 0; // not really required
   }

   for(unsigned int icell=0; icell<n_cell; ++icell)
   {
      auto cell = get_cell_vertices(icell);
      for(unsigned int ipoint=0; ipoint<cell.first; ++ipoint)
      {
         unsigned int gnode = cell.second[ipoint]; // Get the global node number
         unsigned int istor = esup2[gnode]; // location in esup1 where the element will be stored
         esup1[istor] = icell;
         ++esup2[gnode];
      }
   }

   for(unsigned int i=n_vertex; i>0; --i)
   {
      esup2[i] = esup2[i-1]; // reshuffle the esup2
   }
   esup2[0] = 0;
   has_esup = true;
   cout << "Done\n";
}

// Find points surrounding a point
// If all_points==false, find only points connected by an edge.
void Grid::construct_psup(const bool all_points)
{
   if(has_psup) return;
   construct_esup(); // psup needs esup data
   cout << "Constructing points surrounding point ... ";

   unsigned int *lpoin; // Help array to avoid duplication from neighbouring cells
   std::vector<unsigned int> psup1_temp(0);  // temporary vector to allow resize function
   unsigned int istor = 0, gnode=0;
   psup2 = new unsigned int[n_vertex+1];
   lpoin = new unsigned int[n_vertex];

   // Initialize
   for(unsigned int i=0; i<n_vertex; ++i)
   {
      psup2[i] = 0;
      lpoin[i] = 0;
   }
   psup2[n_vertex] = 0;
   
   for(unsigned int ipoint=0; ipoint<n_vertex; ++ipoint) // Loop over all the nodes
   {
      auto esup = get_esup(ipoint); // get cells surrounding the node
      for(unsigned int icell=0; icell<esup.first; ++icell) // Loop over cells surrounding node
      {
         auto cell = get_cell_vertices(esup.second[icell]);
         for(unsigned int jpoint=0; jpoint<cell.first; ++jpoint) // Loop over nodes of cell
         {
            gnode = cell.second[jpoint];  // global node number
            if(gnode != ipoint && lpoin[gnode] != ipoint+1) // check for duplication
            {
               if(all_points == true) // get all the points surrounding the node
               {
                  psup1_temp.resize(psup1_temp.size()+1);
                  psup1_temp[istor] = gnode;
                  lpoin[gnode] = ipoint+1; // using ipoint+1 as 0 is the initialized value
                  ++istor;
               }
               else
               {
                  int prev_node = jpoint-1, next_node = (jpoint+1) % cell.first;
                  if(prev_node < 0) prev_node = cell.first-1;
                  // Check if connected by edge
                  if(cell.second[prev_node] == ipoint || cell.second[next_node] == ipoint)
                  {
                     psup1_temp.resize(psup1_temp.size()+1);
                     psup1_temp[istor] = gnode;
                     lpoin[gnode] = ipoint+1;
                     ++istor;
                  }
               }
            }
         }
      }
      psup2[ipoint+1] = istor;
   }
   // Copy data from psup1_temp to psup1
   psup1 = new unsigned int[psup2[n_vertex]];
   for(unsigned int i=0; i<psup1_temp.size(); ++i)
      psup1[i] = psup1_temp[i];

   psup1_temp.resize(0);
   delete [] lpoin;
   has_psup = true;
   cout << "Done\n";
}
