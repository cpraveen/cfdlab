/* @f$Id: @ref step_1 "step-1".cc 25311 2012-03-20 18:56:09Z bangerth @f$
 *
 * Copyright (C) 1999, 2000, 2001, 2002, 2003, 2005, 2006, 2007, 2009, 2011, 2012 by the deal.II authors
 *
 * This file is subject to QPL and may not be  distributed
 * without copyright and license information. Please refer
 * to the file deal.II/doc/license.html for the  text  and
 * further information on this license.
 */
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>

#include <fstream>
#include <cmath>

using namespace dealii;

void first_grid ()
{
   Triangulation<2> triangulation;
   GridGenerator::hyper_cube (triangulation);
   triangulation.refine_global (4);
   std::ofstream out ("grid-1.eps");
   GridOut grid_out;
   grid_out.write_eps (triangulation, out);
   std::cout << "=====First grid\n";
   std::cout << "Number of cells        = " << triangulation.n_cells() 
             << std::endl;
   std::cout << "Number of active cells = " << triangulation.n_active_cells()
             << std::endl;
}

void second_grid ()
{
   Triangulation<2> triangulation;
   const Point<2> center (1,0);
   const double inner_radius = 0.5,
   outer_radius = 1.0;
   GridGenerator::hyper_shell (triangulation,
                               center, inner_radius, outer_radius,
                               10);
   const HyperShellBoundary<2> boundary_description(center);
   triangulation.set_boundary (0, boundary_description);
   for (unsigned int step=0; step<5; ++step)
   {
      Triangulation<2>::active_cell_iterator
         cell = triangulation.begin_active(),
         endc = triangulation.end();
      for (; cell!=endc; ++cell)
         for (unsigned int v=0;
              v < GeometryInfo<2>::vertices_per_cell;
              ++v)
         {
            const double distance_from_center
               = center.distance (cell->vertex(v));
            if (std::fabs(distance_from_center - inner_radius) < 1e-10)
            {
               cell->set_refine_flag ();
               break;
            }
         }
      triangulation.execute_coarsening_and_refinement ();
      std::cout << "=====Second grid, Step = " << step << std::endl;
      std::cout << "Number of cells        = " << triangulation.n_cells() 
                << std::endl;
      std::cout << "Number of active cells = " << triangulation.n_active_cells()
                << std::endl;
   }
   std::ofstream out ("grid-2.eps");
   GridOut grid_out;
   grid_out.write_eps (triangulation, out);
   triangulation.set_boundary (0);
}

void third_grid ()
{
   Triangulation<2> triangulation;
   GridGenerator::hyper_cube (triangulation);
   triangulation.refine_global (4);
   
   const Point<2> center (0.5,0.5);
   
   for (unsigned int step=0; step<2; ++step)
   {
      Triangulation<2>::active_cell_iterator
         cell = triangulation.begin_active(),
         endc = triangulation.end();
      for (; cell!=endc; ++cell)
         for (unsigned int v=0;
              v < GeometryInfo<2>::vertices_per_cell;
              ++v)
         {
            const double distance_from_center
               = center.distance (cell->vertex(v));
            if (distance_from_center > 0.2 && distance_from_center < 0.4)
            {
               cell->set_refine_flag ();
               break;
            }
         }
      triangulation.execute_coarsening_and_refinement ();
   }
   
   std::ofstream out ("grid-3.eps");
   GridOut grid_out;
   grid_out.write_eps (triangulation, out);
   
   std::ofstream gnu ("grid-3.gpl");
   grid_out.write_gnuplot (triangulation, gnu);
}

int main ()
{
   first_grid ();
   second_grid ();
   third_grid ();
}

