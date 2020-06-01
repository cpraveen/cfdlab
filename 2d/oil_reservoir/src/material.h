#ifndef __MATERIAL_H__
#define __MATERIAL_H__

#include <vector>
#include <cmath>

#define SOLID  0 // no penetration boundary
#define INLET  1 // inlet
#define OUTLET 2 // outlet

#define harmonic_average(a,b)   (2.0*(a)*(b)/((a)+(b)))

extern double viscosity_oil;
extern double density_water;
extern double density_oil;
extern double gravity;

namespace Permeability
{
   extern unsigned int N;
   extern std::vector<double> xl, yl;
}

// viscosity of water
double viscosity_water (const double& concentration);

// mobility of water
double mobility_water (const double& saturation, const double& concentration);

// mobility of oil
double mobility_oil (const double& saturation, const double& concentration);

// total mobility
double mobility_total (const double& saturation, const double& concentration);

// permeability of rock
double rock_permeability (const double& x, const double& y);


// viscosity of water as function of polymer concentration
inline
double viscosity_water (const double& concentration)
{
   return (1.0 + concentration);
}

// mobility of water
inline
double mobility_water (const double& saturation, const double& concentration)
{
   return saturation * saturation / viscosity_water (concentration);
}

// mobility of oil
inline
double mobility_oil (const double& saturation, const double& concentration)
{
   return (1.0 - saturation) * (1.0 - saturation) / viscosity_oil;
}

// total mobility
inline
double mobility_total (const double& saturation, const double& concentration)
{
   return mobility_water (saturation, concentration) +
          mobility_oil   (saturation, concentration);
}

#endif
