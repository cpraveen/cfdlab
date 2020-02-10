#include <iostream>
#include <fstream>
#include <vector>
#include <utility> // std::pair, std::make_pair

using namespace std;

class Grid
{
public:
   Grid ();
   ~Grid ();
   void read_gmsh(const string grid_file);
   void write_vtk(const string grid_file);
   void construct_esup();
   void construct_psup(bool all_points=true);
   void compute_carea();

   unsigned int get_n_vertex()
   {
      return n_vertex;
   }

   unsigned int get_n_cell()
   {
      return n_cell;
   }

   unsigned int get_n_tri()
   {
      return n_tri;
   }

   unsigned int get_n_quad()
   {
      return n_quad;
   }

   unsigned int get_n_bface()
   {
      return n_bface;
   }

   double* get_coord(unsigned int i)
   {
      return &coord[i*dim];
   }

   double get_carea(unsigned int i)
   {
      return carea[i];
   }

   std::pair<unsigned int,const unsigned int*> get_cell_vertices(unsigned int i)
   {
      unsigned int start = cell2[i];
      unsigned int end = cell2[i+1];
      return std::make_pair(end-start,&cell1[start]);
   }

   std::pair<unsigned int,const unsigned int*> get_esup(unsigned int i)
   {
      unsigned int start = esup2[i];
      unsigned int end = esup2[i+1];
      return std::make_pair(end-start,&esup1[start]);
   }

   std::pair<unsigned int,const unsigned int*> get_psup(unsigned int i)
   {
      unsigned int start = psup2[i];
      unsigned int end = psup2[i+1];
      return std::make_pair(end-start,&psup1[start]);
   }

private:
   int          dim;
   unsigned int n_vertex, n_cell, n_tri, n_quad, n_bface;
   double       *coord;
   unsigned int *cell1, *cell2;
   unsigned int *esup1, *esup2;
   unsigned int *psup1, *psup2;
   unsigned int *bface;
   int          *bface_type;
   int          *ctype;
   double       *carea;
   double       *flen;  // length of face
   unsigned int *face;  // face
   double       *fnorm; // unit normal to face
   unsigned int *fnbr;
};
