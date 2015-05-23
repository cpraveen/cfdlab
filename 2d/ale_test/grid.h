#ifndef __GRID_H__
#define __GRID_H__

#include <vector>
#include "vec.h"
#include "face.h"

class Cell
{
   public:
      Vector       centroid;
      unsigned int vertex[3];
};

struct Vertex
{
   Vector coord;
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
      void read_gmsh (std::string grid_file);
      void move_lagrange(double dt);
      void move_noshear(double dt);
      void save();
};

#endif
