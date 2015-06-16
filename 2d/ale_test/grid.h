#ifndef __GRID_H__
#define __GRID_H__

#include <vector>
//#include "parameter.h"
#include "vec.h"
#include "face.h"

class Cell
{
   public:
      Vector       centroid;
      unsigned int vertex[3];
      int          face[3];
      double       area;
      double       radius;
      std::vector<Vector> normal; // Inward normal
};

struct Vertex
{
   Vector coord;
   std::vector<unsigned int> nbr_vertex;
   std::vector<unsigned int> face;
   double radius;
};

class Grid
{
   public:
      Grid () { n_vertex = n_cell = n_face = n_boundary_face = 0; };
      unsigned int n_vertex;
      unsigned int n_cell;
      unsigned int n_face;
      unsigned int n_boundary_face;
      std::vector<Vertex> vertex;
      std::vector<Cell>   cell;
      std::vector<Face>   face;
      std::vector<Face>   bface;
      std::vector<double> mcarea;
      std::vector<double> dcarea;
      std::vector<unsigned int> old_num;
      std::vector<unsigned int> new_num;
      void read_gmsh (std::string grid_file);
      void preproc ();
      void move_lagrange(double dt);
      void move_noshear(double dt);
      void save();

   private:
      void read_bamg (std::string grid_file);
      void read_delaundo (std::string grid_file);
      void compute_cell_centroid ();
      void compute_face_centroid ();
      void compute_cell_area ();
      void compute_face_normal_and_area ();
      void add_face (const Face& new_face);
      void find_vertex_opposite_face ();
      void make_faces ();
      void weight_average () ;
      void vertex_weight_check () ;
      void find_cell_faces ();
      void info ();
      void remove_empty_faces ();
      void renumber();
      void find_nbr_vertex ();
      void print_cells();
      void compute_radius ();

      std::vector< std::vector<unsigned int> > node_face;

};

#endif
