#include <iostream>
#include <cassert>
#include <vector>
#include "grid.h"

using namespace std;

void Grid::allocate ()
{
   assert (nx > 1);
   assert (ny > 1);

   x.allocate (nx+2,ny+2);
   y.allocate (nx+2,ny+2);

   xc.allocate (nx+1,ny+1);
   yc.allocate (nx+1,ny+1);

   ibeg.resize (n_boundary);
   iend.resize (n_boundary);
   jbeg.resize (n_boundary);
   jend.resize (n_boundary);
   boundary_condition.resize (n_boundary);
   b_type.resize (n_boundary);

   n_cells = (nx+1) * (ny+1);
}

// map from (i,j) to cell number
// cells are numbered starting at bottom left corner, then going
// horizontally from i=0 to i=grid.nx+1
unsigned int Grid::cell_num (const unsigned int i, const unsigned int j)
{
   return i + (nx+1)*j;
}
