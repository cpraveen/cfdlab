#include <iostream>
#include <fstream>
#include <sstream>
#include "grid.h"

using namespace std;

Vector gresho(const Vector p)
{
   double r = p.norm();
   double vtheta;
   if(r < 0.2)
      vtheta = 5*r;
   else if(r > 0.4)
      vtheta = 0;
   else
      vtheta = 2 - 5*r;

   double theta = atan2(p.y, p.x);
   Vector vel;
   vel.x = -vtheta * sin(theta);
   vel.y =  vtheta * cos(theta);
   return vel;
}

void gresho_gradu(const Vector& p, double gradu[2][2])
{
   double r = p.norm();
   if(r < 0.2)
   {
      gradu[0][0] =  0.0; gradu[0][1] = -5.0;
      gradu[1][0] =  5.0; gradu[1][1] =  0.0;
   }
   else if(r > 0.4)
   {
      gradu[0][0] =  0.0; gradu[0][1] =  0.0;
      gradu[1][0] =  0.0; gradu[1][1] =  0.0;
   }
   else
   {
      double r3 = r*r*r;
      gradu[0][0] =  2.0*p.x*p.y/r3        ; gradu[0][1] = -2/r + 2*p.y*p.y/r3 + 5;
      gradu[1][0] =  2/r - 2*p.x*p.x/r3 + 5; gradu[1][1] = -2*p.x*p.y/r3;
   }

}

void Grid::move_lagrange(double dt)
{
   for(unsigned int i=0; i<n_vertex; ++i)
   {
      Vector vel = gresho(vertex[i].coord);
      vertex[i].coord += vel * dt;
   }
}

void Grid::move_noshear(double dt)
{
   vector<Vector> vel (n_vertex);
   vector<int> nbr (n_vertex);
   for(unsigned int i=0; i<n_vertex; ++i)
   {
      vel[i] = 0;
      nbr[i] = 0;
   }

   double gradu[2][2], S[2][2], A[2][2];
   for(unsigned int i=0; i<n_cell; ++i)
   {
      gresho_gradu(cell[i].centroid, gradu);
      for(unsigned int m=0; m<2; ++m)
         for(unsigned int n=0; n<2; ++n)
         {
            S[m][n] = 0.5*(gradu[m][n] + gradu[n][m]);
            A[m][n] = 0.5*(gradu[m][n] - gradu[n][m]);
         }
      double div = S[0][0] + S[1][1];
      A[0][0] += 0.5 * div;
      A[1][1] += 0.5 * div;

      for(unsigned int j=0; j<3; ++j)
      {
         int n = cell[i].vertex[j];
         Vector dr = vertex[n].coord - cell[i].centroid;
         Vector velc = gresho(cell[i].centroid);
         velc.x += A[0][0] * dr.x + A[0][1] * dr.y;
         velc.y += A[1][0] * dr.x + A[1][1] * dr.y;
         vel[n] += velc;
         nbr[n] += 1;
      }
   }

   for(unsigned int i=0; i<n_vertex; ++i)
   {
      vel[i] *= (1.0/nbr[i]);
      vertex[i].coord += vel[i] * dt;
   }
}

void Grid::save()
{
   static int counter = 0;

   string filename;
   if     (counter <= 9)    filename = "sol000";
   else if(counter <= 99)   filename = "sol00";
   else if(counter <= 999)  filename = "sol0";
   else if(counter <= 9999) filename = "sol";
   else
   {
      cout << "Writer::output: counter is too large !!!\n";
      exit(0);
   }
   stringstream ss;
   ss << counter;
   filename += ss.str() + ".vtk";

   ofstream vtk;
   vtk.open (filename.c_str());

   vtk << "# vtk DataFile Version 3.0" << endl;
   vtk << "flo3d" << endl;
   vtk << "ASCII" << endl;
   vtk << "DATASET UNSTRUCTURED_GRID" << endl;
   vtk << "POINTS  " << n_vertex << "  float" << endl;

   for(unsigned int i=0; i<n_vertex; ++i)
      vtk << vertex[i].coord.x << " " 
          << vertex[i].coord.y << " " 
          << 0.0 << endl;

   vtk << "CELLS  " << n_cell << " " << 4 * n_cell << endl;
   for(unsigned int i=0; i<n_cell; ++i)
      vtk << 3 << "  " 
          << cell[i].vertex[0] << " "
          << cell[i].vertex[1] << " "
          << cell[i].vertex[2] << endl;

   vtk << "CELL_TYPES  " << n_cell << endl;
   for(unsigned int i=0; i<n_cell; ++i)
      vtk << 5 << endl;

   vtk.close();

   ++counter;
}

int main(int argc, char* argv[])
{
   Grid grid;
   grid.read_gmsh("square.msh");
   grid.save();

   double dt = 0.01;
   for(unsigned int i=0; i<20; ++i)
   {
      //grid.move_lagrange(dt);
      grid.move_noshear(dt);
      grid.save();
   }
}
