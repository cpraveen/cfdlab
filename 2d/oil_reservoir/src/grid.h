#ifndef __GRID_H__
#define __GRID_H__

#include <vector>
#include "matrix.h"

enum BType { imin, imax, jmin, jmax };

// class to hold the grid object
class Grid
{
   public:

      Grid () { nx = ny = n_cells = 0; };
      ~Grid () {};

      unsigned int nx, ny;
      Matrix x, y;   // grid vertices
      Matrix xc, yc; // cell centers
      unsigned int n_boundary;
      unsigned int n_cells;
      double xmin, xmax, ymin, ymax;
      double dx, dy;

      unsigned int cell_num(const unsigned int, const unsigned int);

      std::vector<unsigned int> ibeg;
      std::vector<unsigned int> iend;
      std::vector<unsigned int> jbeg;
      std::vector<unsigned int> jend;
      std::vector<unsigned int> boundary_condition;
      std::vector<BType> b_type;

      void allocate ();

};

#endif
