#include <iostream>
#include <cmath>
#include <cassert>
#include <fstream>
#include <cstdlib>
#include "grid.h"

extern bool debug;

using namespace std;

//------------------------------------------------------------------------------
// Compute cell Centroid
//------------------------------------------------------------------------------
void Grid::compute_cell_centroid ()
{

      cout << " median cell\n";
      for(unsigned int i=0; i<n_cell; ++i)
      { 
         unsigned int v0, v1, v2;
         v0 = cell[i].vertex[0];
         v1 = cell[i].vertex[1];
         v2 = cell[i].vertex[2];
         cell[i].centroid = ( vertex[v0].coord + 
                              vertex[v1].coord + 
                              vertex[v2].coord ) / 3.0;
      }   
}

//------------------------------------------------------------------------------
// Compute face Centroid
//------------------------------------------------------------------------------
void Grid::compute_face_centroid ()
{
   cout << "Computing face centers ...\n";
   for(unsigned int i=0; i<n_face; ++i)
   { 
      unsigned int v0 = face[i].vertex[0];
      unsigned int v1 = face[i].vertex[1];
      face[i].centroid = ( vertex[v0].coord + vertex[v1].coord ) / 2.0;
   }   

   for(unsigned int i=0; i<bface.size(); ++i)
   { 
      unsigned int v0 = bface[i].vertex[0];
      unsigned int v1 = bface[i].vertex[1];
      bface[i].centroid = ( vertex[v0].coord + vertex[v1].coord ) / 2.0;
   }   
}
//------------------------------------------------------------------------------
// Compute cell areas
//------------------------------------------------------------------------------
void Grid::compute_cell_area ()
{
   cout << "Computing cell areas ...\n";

   mcarea.resize (n_vertex);
   dcarea.resize (n_vertex);

   for(unsigned int i=0; i<n_vertex; ++i)
   {
      mcarea[i] = 0.0;
      dcarea[i] = 0.0;
   }

   for(unsigned int i=0; i<n_cell; ++i)
   {
      unsigned int v0 = cell[i].vertex[0];
      unsigned int v1 = cell[i].vertex[1];
      unsigned int v2 = cell[i].vertex[2];

      // compute area as vector cross product
      cell[i].area = (vertex[v1].coord - vertex[v0].coord) ^ 
                     (vertex[v2].coord - vertex[v0].coord);
      cell[i].area *= 0.5;

      assert ( cell[i].area > 0.0 );

      // contribution to dual cell
      mcarea[v0] += cell[i].area / 3.0;
      mcarea[v1] += cell[i].area / 3.0;
      mcarea[v2] += cell[i].area / 3.0;
   }

      dcarea = mcarea;

}

//------------------------------------------------------------------------------
// Compute face normals
//------------------------------------------------------------------------------
void Grid::compute_face_normal_and_area ()
{

   // Check orientation of interior faces
   for(unsigned int i=0; i<n_face; ++i)
   {
      unsigned int v0 = face[i].vertex[0];
      unsigned int v1 = face[i].vertex[1];
      Vector dr = vertex[v1].coord - vertex[v0].coord;
      
      Vector normal;
      normal.x = +dr.y;
      normal.y = -dr.x;
            
      // Check orintation of boundary face
      unsigned int cl = face[i].lcell;
      Vector dcf = face[i].centroid - cell[cl].centroid;
      double d = dcf * normal;
      if(d < 0.0)
      {
         face[i].vertex[0] = v1;
         face[i].vertex[1] = v0;
      }
   }
   
   // interior faces
   for(unsigned int i=0; i<n_face; ++i)
   {
      unsigned int cl = face[i].lcell;
      Vector dr = face[i].centroid - cell[cl].centroid;

      if(face[i].type == -1) // interior edge, has right cell also
      {
         unsigned int cr = face[i].rcell;
         dr += cell[cr].centroid - face[i].centroid;
      }

      face[i].normal.x = -dr.y;
      face[i].normal.y =  dr.x;

      face[i].measure = face[i].normal.norm();
  }

   // boundary faces
   for(unsigned int i=0; i<bface.size(); ++i)
   {
      unsigned int v0 = bface[i].vertex[0];
      unsigned int v1 = bface[i].vertex[1];
      Vector dr = vertex[v1].coord - vertex[v0].coord;

      bface[i].normal.x =  dr.y;
      bface[i].normal.y = -dr.x;

      bface[i].measure = bface[i].normal.norm();
      
      // Check orintation of boundary face
      unsigned int cl = bface[i].lcell;
      Vector dcf = bface[i].centroid - cell[cl].centroid;
      double d = dcf * bface[i].normal;
      if(d < 0.0)
      {
         bface[i].vertex[0] = v1;
         bface[i].vertex[1] = v0;
         bface[i].normal *= -1.0;
      }
   }

   // Gradient of P1 basis function on each triangle
   // You need to divide by 2*area to get gradient
   Vector dr01, dr12, dr20;
   for(unsigned int i=0; i<n_cell; ++i)
   {
      cell[i].normal.resize(3);

      unsigned int v0 = cell[i].vertex[0];
      unsigned int v1 = cell[i].vertex[1];
      unsigned int v2 = cell[i].vertex[2];

      dr01.x =  (vertex[v0].coord.y - vertex[v1].coord.y);
      dr12.x =  (vertex[v1].coord.y - vertex[v2].coord.y);
      dr20.x =  (vertex[v2].coord.y - vertex[v0].coord.y);

      dr01.y = -(vertex[v0].coord.x - vertex[v1].coord.x);
      dr12.y = -(vertex[v1].coord.x - vertex[v2].coord.x);
      dr20.y = -(vertex[v2].coord.x - vertex[v0].coord.x);

      cell[i].normal[0] = dr12;
      cell[i].normal[1] = dr20;
      cell[i].normal[2] = dr01;
   }


}

//------------------------------------------------------------------------------
// Add new_face to face list
//------------------------------------------------------------------------------
void Grid::add_face (const Face& new_face)
{
   bool found = false;
   unsigned int n = 0;

   // Take any vertex of this face
   unsigned int v = new_face.vertex[0];

   // Loop over all existing faces of vertex v
   while (n<node_face[v].size() && !found)
   {
      unsigned int f = node_face[v][n];

      if(face[f].lcell==-1) // Boundary face not filled yet
      {
         if(face[f] == new_face)
         {
            face[f].lcell   = new_face.lcell;
            found = true;
         }
      }
      else if(face[f].rcell == -1) // Boundary or interior face
      {
         if(face[f] == new_face)
         {
            face[f].rcell   = new_face.lcell;
            found = true;
         }
      }

      ++n;
   }

   if(!found) // This is a new face
   {
      face.resize (n_face+1);
      face[n_face].type = -1; // TODO: NEED TO GIVE NEW TYPE FOR INTERIOR FACES
      face[n_face].vertex [0] = new_face.vertex [0];
      face[n_face].vertex [1] = new_face.vertex [1];
      face[n_face].lcell      = new_face.lcell; // left cell
      face[n_face].rcell      = -1; // right cell to be found

      // Add this face to its two vertices
      for(unsigned int i=0; i<2; ++i)
      {
         v = new_face.vertex[i];
         node_face[v].push_back (n_face);
      }

      ++n_face;
   }

}

//------------------------------------------------------------------------------
// Find vertex opposite each face; for boundary faces only left vertex present
//------------------------------------------------------------------------------
void Grid::find_vertex_opposite_face ()
{
   for(unsigned int i=0; i<n_face; ++i)
      {
         int vl = -1;
         int v0 = face[i].vertex[0];
         int v1 = face[i].vertex[1];
         int cl = face[i].lcell;
         for(unsigned int j=0; j<3; ++j)
            if(cell[cl].vertex[j] != v0 && cell[cl].vertex[j] != v1)
               vl = cell[cl].vertex[j];
         assert (vl != -1);
         face[i].lvertex = vl;

         // interior face: find right vertex also
         if(face[i].type == -1)
         {
            int vr = -1;
            int cr = face[i].rcell;
            for(unsigned int j=0; j<3; ++j)
               if(cell[cr].vertex[j] != v0 && cell[cr].vertex[j] != v1)
                  vr = cell[cr].vertex[j];
            assert (vr != -1);
            face[i].rvertex = vr;
         }
         else
            face[i].rvertex = -1;
      }
}

//------------------------------------------------------------------------------
// Create interior faces and connectivity data
//------------------------------------------------------------------------------
void Grid::make_faces ()
{
   cout << "Creating faces ..." << endl;
   unsigned int i;

   node_face.resize (n_vertex);

   // Existing boundary faces
   for(i=0; i<n_face; ++i)
   {
      face[i].lcell   = -1;
      face[i].rcell   = -1;

      // Add this face to the two vertices
      for(unsigned int j=0; j<2; ++j)
      {
         unsigned int v = face[i].vertex[j];
         node_face[v].push_back(i);
      }
   }

   Face new_face;

   for(i=0; i<n_cell; ++i)
   {
      // first face
      new_face.vertex[0] = cell[i].vertex[0];
      new_face.vertex[1] = cell[i].vertex[1];
      new_face.lcell     = i;
      add_face (new_face);

      // second face
      new_face.vertex[0] = cell[i].vertex[1];
      new_face.vertex[1] = cell[i].vertex[2];
      new_face.lcell     = i;
      add_face (new_face);

      // third face
      new_face.vertex[0] = cell[i].vertex[2];
      new_face.vertex[1] = cell[i].vertex[0];
      new_face.lcell     = i;
      add_face (new_face);
   }

   cout << "Checking face data ..." << endl;

   // Now check that face data is complete
   for(i=0; i<n_face; ++i)
   {
      // Check left cell exists
      assert(face[i].lcell != -1);

      if(face[i].type == -1) // Interior face, check right cell
         assert(face[i].rcell != -1);
   }

   find_vertex_opposite_face ();

   // Copy boundary faces into bface
   for(i=0; i<n_face; ++i)
      if(face[i].type != -1)
         bface.push_back(face[i]);

   // Free memory of node_face since we dont need it any more
   for(i=0; i<n_vertex; ++i)
      node_face[i].resize (0);
   node_face.resize (0);
}

//------------------------------------------------------------------------------
// Find cells surrounding a cell
//------------------------------------------------------------------------------
void Grid::find_cell_faces ()
{
   cout << "Finding faces for each cell ..." << endl;
   
   unsigned int i, j;
   int lcell, rcell;
   
   // First put all faces to -1
   for(i=0; i<n_cell; ++i)
   {
      cell[i].face[0] = -1;
      cell[i].face[1] = -1;
      cell[i].face[2] = -1;
   }
   
   for(i=0; i<n_face; ++i)
   {
      lcell = face[i].lcell;
      j = 0;
      while(cell[lcell].face[j] != -1)
         ++j;
      cell[lcell].face[j] = i;
            
      rcell = face[i].rcell;
      if(rcell != -1)
      { 
         j = 0;
         while(cell[rcell].face[j] != -1)
            ++j;
         cell[rcell].face[j] = i;
      }
    }
   
}

//------------------------------------------------------------------------------
// Find points connected to a point
//------------------------------------------------------------------------------
void Grid::remove_empty_faces()
{
   cout << "Removing faces with zero area ...\n";
   vector<Face> tmp;
   for(unsigned int i=0; i<n_face; ++i)
   {
      if(face[i].measure > 0.0) tmp.push_back (face[i]);
   }
   cout << "   faces with zero area = " << face.size() - tmp.size() << "\n";

   face.resize (0);
   face.resize (tmp.size());
   face = tmp;
   tmp.resize(0);

   n_face = face.size ();
}

//------------------------------------------------------------------------------
// Find points connected to a point
//------------------------------------------------------------------------------
void Grid::find_nbr_vertex()
{
   cout << "Finding vertices around a vertex ...\n";
   for(unsigned int i=0; i<n_face; ++i)
   {
      unsigned int v0 = face[i].vertex[0];
      unsigned int v1 = face[i].vertex[1];
      vertex[v0].nbr_vertex.push_back (v1);
      vertex[v1].nbr_vertex.push_back (v0);
      vertex[v0].face.push_back (i);
      vertex[v1].face.push_back (i);
   }
}

//------------------------------------------------------------------------------
// Renumbering according to cuthill-mckee algorithm
//------------------------------------------------------------------------------
void Grid::renumber()
{  
   cout << "Renumbering vertices using Cuthill-McKee algorithm ...\n";

   unsigned int i, j, current_cell, old_cell, k;
   int neighbour =-1;
   // new_num vector says directly the value of renumbering tag for a old cell 
   // number old_num vector says what is the value of old cell number for a 
   // given renumbering tag renumbering cell vector is a dummy vector to 
   // reshuffle all old cell according to new numbering 
   new_num.resize(n_vertex,0);
   old_num.resize(n_vertex,0);
   
   // Write initial cell numbering to file
   if(debug)
   {
      ofstream out("number_old.dat");
      
      for(i=0; i<n_vertex; ++i)
      { 
         out << i << " " << i << endl;
         for(j=0; j<vertex[i].nbr_vertex.size(); ++j)
         {    
            neighbour = vertex[i].nbr_vertex[j];
            out << i << " " << neighbour << endl;
         }
      }
      out.close();
   }
   
   k=1; // k is the renumbering tag according to the algorithm
   for(i=0; i<n_vertex; ++i)
   { 
      j = 0;
      current_cell = old_num[i] ;
      while(j < vertex[current_cell].nbr_vertex.size())
      {     
         neighbour = vertex[current_cell].nbr_vertex[j];
         old_cell = neighbour;
         if (new_num[old_cell] == 0 && old_cell != 0)
         {
            old_num[k] = old_cell;
            new_num[old_cell] = k;
            ++k;
         }
         ++j;
      }
   }
   
   // Save new cell numbering to file
   if(debug)
   {
      ofstream out("number_new.dat");
      for(i=0; i<n_vertex; ++i)
      {
         current_cell = old_num[i];
         out << i << " " << i << endl;
         for(j=0; j<vertex[current_cell].nbr_vertex.size(); ++j)
         {    
            neighbour = vertex[current_cell].nbr_vertex[j];
            out << i << " " << new_num[neighbour] << endl;
         }
         
      }
      out.close();
   }
}   

//------------------------------------------------------------------------------
// Print triangles and dual cells into file
// Can be visualized with gnuplot
//------------------------------------------------------------------------------
void Grid::print_cells ()
{
   cout << "Writing triangles into tri.dat ...\n";
   cout << "Writing dual cells into dual.dat ...\n";

   ofstream tri("tri.dat");
   ofstream dual("dual.dat");

   for(unsigned int i=0; i<n_cell; ++i)
   {
      unsigned int n0 = cell[i].vertex[0];
      unsigned int n1 = cell[i].vertex[1];
      unsigned int n2 = cell[i].vertex[2];
      tri << vertex[n0].coord.x << "  " << vertex[n0].coord.y << "\n";
      tri << vertex[n1].coord.x << "  " << vertex[n1].coord.y << "\n";
      tri << vertex[n2].coord.x << "  " << vertex[n2].coord.y << "\n";
      tri << vertex[n0].coord.x << "  " << vertex[n0].coord.y << "\n";
      tri << "\n";

      Vector p01 = (vertex[n0].coord + vertex[n1].coord)/2.0;
      Vector p12 = (vertex[n1].coord + vertex[n2].coord)/2.0;
      Vector p20 = (vertex[n2].coord + vertex[n0].coord)/2.0;

      dual << p01.x << "  " << p01.y << "\n";
      dual << cell[i].centroid.x << "  " << cell[i].centroid.y << "\n";
      dual << "\n";

      dual << p12.x << "  " << p12.y << "\n";
      dual << cell[i].centroid.x << "  " << cell[i].centroid.y << "\n";
      dual << "\n";

      dual << p20.x << "  " << p20.y << "\n";
      dual << cell[i].centroid.x << "  " << cell[i].centroid.y << "\n";
      dual << "\n";
   }
   for(unsigned int i=0; i<bface.size(); ++i)
   {
      unsigned int v0 = bface[i].vertex[0];
      unsigned int v1 = bface[i].vertex[1];
      dual << vertex[v0].coord.x << "  " << vertex[v0].coord.y << endl;
      dual << vertex[v1].coord.x << "  " << vertex[v1].coord.y << endl;
      dual << endl;
   }
   tri.close();
   dual.close();
}

//------------------------------------------------------------------------------
// Preprocess the grid
//------------------------------------------------------------------------------
void Grid::preproc ()
{
   make_faces ();
   find_cell_faces ();
   compute_face_centroid ();
   compute_cell_centroid ();
   compute_cell_area ();
   compute_face_normal_and_area ();
   find_nbr_vertex ();
   renumber ();
   if(debug)
      print_cells();
}
