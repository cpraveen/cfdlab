#include <iostream>
#include <cmath>
#include "material.h"


// permeability of rock
double rock_permeability (const double& x, const double& y)
{
   //return 1.0;
   //return 1.0 + 0.5 * cos(6.0*M_PI*(x+0.2)) * cos(6.0*M_PI*y);

   // random permeability
   double k = 0.0;
   for(unsigned int i=0; i<Permeability::N; ++i)
   {
      k += exp(-(pow(x-Permeability::xl[i],2) + pow(y-Permeability::yl[i],2))/pow(0.05,2));
   }
   k = std::max(k, 0.5);
   k = std::min(k, 1.5);
   return k;
}
