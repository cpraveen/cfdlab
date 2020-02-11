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
   void construct_psup(const bool all_points=true);
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

   const double* get_coord(unsigned int i)
   {
      return &coord[i*dim];
   }

   double get_cell_area(unsigned int i)
   {
      return carea[i];
   }

   const double* get_face_normal(unsigned int i)
   {
      return &fnorm[i*dim];
   }

   const double* get_bface_normal(unsigned int i)
   {
      return &bface_norm[i*dim];
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
   const int    dim = 2;
   unsigned int n_vertex, n_cell, n_tri, n_quad, n_bface;
   double       *coord;

   // cell data
   unsigned int *cell1, *cell2;
   int          *ctype;
   double       *carea;

   // connectivity information
   unsigned int *esup1, *esup2;
   unsigned int *psup1, *psup2;

   // boundary face data
   unsigned int *bface;      // vertices forming the face
   unsigned int *bface_cell; // cell adjacent to boundary face
   int          *bface_type; // type read from grid file, used for bc
   double       *bface_norm; // unit outward normal

   // Interior faces
   double       *flen;  // length of face
   unsigned int *face;  // vertex numbers for each face
   double       *fnorm; // unit normal to face
   unsigned int *fcell; // cells adjacent to face

   bool         has_esup;
   bool         has_psup;
};
