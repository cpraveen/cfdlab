#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

void write_rectilinear_grid(int nx,
                            int ny,
                            double *x,
                            double *y,
                            double **var)
{
   int nz = 1; // We have a 2d grid

   ofstream fout;
   fout.open("rect.vtk");
   fout << "# vtk DataFile Version 2.0" << endl;
   fout << "Cartesian grid" << endl;
   fout << "ASCII" << endl;
   fout << "DATASET RECTILINEAR_GRID" << endl;
   fout << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
   fout << "X_COORDINATES " << nx << " float" << endl;
   for(int i=0; i<nx; ++i)
      fout << x[i] << " ";
   fout << endl;
   fout << "Y_COORDINATES " << ny << " float" << endl;
   for(int j=0; j<ny; ++j)
      fout << y[j] << " ";
   fout << endl;
   fout << "Z_COORDINATES " << nz << " float" << endl;
   fout << 0.0 << endl;

   fout << "POINT_DATA " << nx*ny*nz << endl;
   fout << "SCALARS density float" << endl;
   fout << "LOOKUP_TABLE default" << endl;
   for(int i=0; i<nx; ++i)
   {
      for(int j=0; j<ny; ++j)
         fout << var[i][j] << " ";
      fout << endl;
   }
   fout.close();

   cout << "Wrote Cartesian grid into rect.vtk" << endl;
}

void test_rectilinear_grid()
{
   cout << "Saving rectilinear grid\n";
   const int nx = 100, ny = 50;
   const double xmin = 0.0, xmax = 2.0;
   const double ymin = 0.0, ymax = 1.0;
   const double dx = (xmax-xmin)/(nx-1);
   const double dy = (ymax-ymin)/(ny-1);
   double *x, *y, **var;

   x = new double[nx];
   y = new double[ny];

   for(int i=0; i<nx; ++i)
      x[i] = xmin + i*dx;

   for(int j=0; j<ny; ++j)
      y[j] = ymin + j*dy;

   var = new double*[nx];
   for(int i=0; i<nx; ++i)
      var[i] = new double[ny];

   for(int i=0; i<nx; ++i)
      for(int j=0; j<ny; ++j)
      {
         var[i][j] = sin(2.0*M_PI*x[i]) * sin(2.0*M_PI*y[j]);
      }
   write_rectilinear_grid(nx, ny, x, y, var);

   // Deallocate memory here
}

void write_structured_grid(int ni,
                           int nj,
                           double **x,
                           double **y,
                           double **var)
{
   int nk = 1; // We have a 2d grid

   ofstream fout;
   fout.open("struct.vtk");
   fout << "# vtk DataFile Version 2.0" << endl;
   fout << "Structured grid" << endl;
   fout << "ASCII" << endl;
   fout << "DATASET STRUCTURED_GRID" << endl;
   fout << "DIMENSIONS " << ni << " " << nj << " " << nk << endl;
   fout << "POINTS " << ni*nj*nk << " float" << endl;
   for(int i=0; i<ni; ++i)
      for(int j=0; j<nj; ++j)
         fout << x[i][j] << " " << y[i][j] << " 0.0" << endl;

   fout << "POINT_DATA " << ni*nj*nk << endl;
   fout << "SCALARS density float" << endl;
   fout << "LOOKUP_TABLE default" << endl;
   for(int i=0; i<ni; ++i)
   {
      for(int j=0; j<nj; ++j)
         fout << var[i][j] << " ";
      fout << endl;
   }
   fout.close();

   cout << "Wrote structured grid into struct.vtk" << endl;
}

void test_structured_grid()
{
   cout << "Saving structured grid\n";
   const int nt = 100, nr = 100;
   const double rmin = 1.0, rmax = 2.0;
   const double tmin = 0.0, tmax = 0.5*M_PI;
   const double dr = (rmax - rmin)/(nr-1);
   const double dt = (tmax - tmin)/(nt-1);

   double **x, **y, **var;

   x = new double*[nt];
   y = new double*[nt];
   var = new double*[nt];
   for(int i=0; i<nt; ++i)
   {
      x[i] = new double[nr];
      y[i] = new double[nr];
      var[i] = new double[nr];
   }

   for(int i=0; i<nt; ++i)
      for(int j=0; j<nr; ++j)
      {
         double t = tmin + i*dt;
         double r = rmin + j*dr;
         x[i][j] = r * sin(t);
         y[i][j] = r * cos(t);
         var[i][j] = sin(2.0*M_PI*x[i][j]) * sin(2.0*M_PI*y[i][j]);
      }

   write_structured_grid(nt, nr, x, y, var);
}

int main()
{
   test_rectilinear_grid();
   test_structured_grid();
}
