#include <iostream>
#include "material.h"
#include "reservoir.h"

// Global variables
double gravity;
double pinlet;
double poutlet;
double viscosity_oil;
double density_water;
double density_oil;

unsigned int Permeability::N;
std::vector<double> Permeability::xl;
std::vector<double> Permeability::yl;

// counter to count how many times the flux has an interior minimum
int n_interior_min;

using namespace std;

int main ()
{
   cout << "Starting reservoir problem ..." << endl;

   ReservoirProblem  reservoir_problem;

   reservoir_problem.run ();

   return 0;
}
