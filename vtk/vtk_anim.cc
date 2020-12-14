#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

using namespace std;

string get_filename(string base_name,
                    int ndigits,
                    int c)
{
   if(c > pow(10,ndigits))
   {
      cout << "get_filename: Not enough digits !!!\n";
      cout << "ndigits= " << ndigits << endl;
      cout << "c      = " << c << endl;
      exit(0);
   }

   string name = base_name;
   // Improve the below part
   if(c < 10)
   {
      for(int i=0; i<ndigits-1; ++i)
         name += "0";
   }
   else if(c < 100)
   {
      for(int i=0; i<ndigits-2; ++i)
         name += "0";
   }
   else if(c < 1000)
   {
      for(int i=0; i<ndigits-3; ++i)
         name += "0";
   }
   else if(c < 10000)
   {
      for(int i=0; i<ndigits-4; ++i)
         name += "0";
   }
   else
   {
      cout << "get_filename: Not implemented\n";
      exit(0);
   }

   name += to_string(c) + ".vtk";
   return name;
}

void write_rectilinear_grid(int nx,
                            int ny,
                            double *x,
                            double *y,
                            double **var,
                            double t,
                            int c,
                            string filename)
{
   int nz = 1; // We have a 2d grid

   ofstream fout;
   fout.open(filename);
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

   cout << "Wrote Cartesian grid into " << filename << endl;
}

int main()
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

   double u = 1.0, v = 1.0;
   double t = 0.0, dt = 0.1, Tf = 10.0;
   int c = 0;
   while(t < Tf)
   {
      for(int i=0; i<nx; ++i)
         for(int j=0; j<ny; ++j)
         {
            double xx = x[i] - u*t;
            double yy = y[j] - v*t;
            var[i][j] = sin(2.0*M_PI*xx) * sin(2.0*M_PI*yy);
         }
      string filename = get_filename("sol_",3,c);
      write_rectilinear_grid(nx, ny, x, y, var, t, c, filename);
      t += dt;
      ++c;
   }

   // Deallocate memory here
   delete[] x;
   delete[] y;
   for(int i=0; i<nx; ++i)
      delete[] var[i];
   delete[] var;

   return 0;
}
