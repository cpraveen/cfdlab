#include <math.h>

#include "libmesh/libmesh_common.h"

using namespace libMesh;

Real exact_solution (const Real x,
                     const Real y,
                     const Real t)
{
   Real r2 = std::pow(x-0.5, 2) + std::pow(y, 2);
   Real r  = std::sqrt(r2);

   if(r < 0.25)
      return std::pow( std::cos(2*M_PI*r), 2 );
   else
      return 0.0;
}
