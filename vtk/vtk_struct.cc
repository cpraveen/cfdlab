#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

void write_rectilinear_grid(int nx,
                            int ny,
                            double *x,
                            double *y,
                            double **var,
                            double t,
                            int c)
{
   int nz = 1; // We have a 2d grid

   ofstream fout;
   fout.open("rect.vtk");
   fout << "# vtk DataFile Version 3.0" << endl;
   fout << "Cartesian grid" << endl;
   fout << "ASCII" << endl;
   fout << "DATASET RECTILINEAR_GRID" << endl;
   fout << "FIELD FieldData 2" << endl;
   fout << "TIME 1 1 double" << endl;
   fout << t << endl;
   fout << "CYCLE 1 1 int" << endl;
   fout << c << endl;
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
   // no need for k-loop since nk=1
   for(int j=0; j<ny; ++j)
   {
      for(int i=0; i<nx; ++i)
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
   write_rectilinear_grid(nx, ny, x, y, var, 0.0, 0);

   // Deallocate memory here
   delete[] x;
   delete[] y;
   for(int i=0; i<nx; ++i)
      delete[] var[i];
   delete[] var;
}

void write_structured_grid(int ni,
                           int nj,
                           double **x,
                           double **y,
                           double **var,
                           double t,
                           int c)
{
   int nk = 1; // We have a 2d grid

   ofstream fout;
   fout.open("struct.vtk");
   fout << "# vtk DataFile Version 3.0" << endl;
   fout << "Structured grid" << endl;
   fout << "ASCII" << endl;
   fout << "DATASET STRUCTURED_GRID" << endl;
   fout << "FIELD FieldData 2" << endl;
   fout << "TIME 1 1 double" << endl;
   fout << t << endl;
   fout << "CYCLE 1 1 int" << endl;
   fout << c << endl;
   fout << "DIMENSIONS " << ni << " " << nj << " " << nk << endl;
   fout << "POINTS " << ni*nj*nk << " float" << endl;
   // no need for k-loop since nk=1
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

void write_tecplot_grid(int ni,
                        int nj,
                        double **x,
                        double **y,
                        double **var,
                        double t,
                        int c)
{
   ofstream fout;
   fout.open("struct.plt");
   fout << "TITLE = \"Tecplot format\"" << endl;
   fout << "VARIABLES = x, y, density" << endl;
   fout << "ZONE STRANDID=1, SOLUTIONTIME=" << t << ", I=" << ni << ", J=" << nj
        << ", DATAPACKING=POINT" << endl;
   for(int j=0; j<nj; ++j)
      for(int i=0; i<ni; ++i)
         fout << x[i][j] << " " << y[i][j] << " " << var[i][j] << endl;

   cout << "Wrote structured grid into struct.plt" << endl;
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
         x[i][j] = r * cos(t);
         y[i][j] = r * sin(t);
         var[i][j] = sin(2.0*M_PI*x[i][j]) * sin(2.0*M_PI*y[i][j]);
      }

   write_structured_grid(nt, nr, x, y, var, 0.0, 0);
   write_tecplot_grid(nt, nr, x, y, var, 0.0, 0);

   // deallocate memory
   for(int i=0; i<nt; ++i)
   {
      delete[] x[i];
      delete[] y[i];
      delete[] var[i];
   }
   delete[] x;
   delete[] y;
   delete[] var;
}

int main()
{
   test_rectilinear_grid();
   test_structured_grid();
   return 0;
}
