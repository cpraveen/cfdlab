/* Author: Wolfgang Bangerth, University of Heidelberg, 1999 */
/*    @f$Id: @ref step_2 "step-2".cc 25560 2012-05-29 19:16:36Z bangerth @f$ */
/*    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2006, 2008, 2009, */
/*                  2010, 2011, 2012 by the deal.II authors       */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>

#include <fstream>

using namespace dealii;

void make_grid (Triangulation<2> &triangulation)
{
   const Point<2> center (1,0);
   const double inner_radius = 0.5,
   outer_radius = 1.0;
   GridGenerator::hyper_shell (triangulation,
                               center, inner_radius, outer_radius,
                               10);
   static const HyperShellBoundary<2> boundary_description(center);
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
   }
}

void distribute_dofs (DoFHandler<2> &dof_handler)
{
   static const FE_Q<2> finite_element(1);
   dof_handler.distribute_dofs (finite_element);
   CompressedSparsityPattern compressed_sparsity_pattern(dof_handler.n_dofs(),
                                                         dof_handler.n_dofs());
   DoFTools::make_sparsity_pattern (dof_handler, compressed_sparsity_pattern);
   SparsityPattern sparsity_pattern;
   sparsity_pattern.copy_from (compressed_sparsity_pattern);
   std::ofstream out ("sparsity_pattern.1");
   sparsity_pattern.print_gnuplot (out);
}

void renumber_dofs (DoFHandler<2> &dof_handler)
{
   DoFRenumbering::Cuthill_McKee (dof_handler);
   CompressedSparsityPattern compressed_sparsity_pattern(dof_handler.n_dofs(),
                                                         dof_handler.n_dofs());
   DoFTools::make_sparsity_pattern (dof_handler, compressed_sparsity_pattern);
   SparsityPattern sparsity_pattern;
   sparsity_pattern.copy_from (compressed_sparsity_pattern);
   std::ofstream out ("sparsity_pattern.2");
   sparsity_pattern.print_gnuplot (out);
}

int main ()
{
   Triangulation<2> triangulation;
   make_grid (triangulation);
   DoFHandler<2> dof_handler (triangulation);
   distribute_dofs (dof_handler);
   renumber_dofs (dof_handler);
}   
