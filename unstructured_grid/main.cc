#include "grid.h"

int main()
{
   Grid grid;
   grid.read_gmsh("cylinder.msh");
   grid.write_vtk("cylinder.vtk");
   grid.construct_esup();
}
