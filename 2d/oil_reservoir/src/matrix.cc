#include <cassert>
#include "matrix.h"

using namespace std;

// empty constructor
Matrix::Matrix ()
{
   nrow = 0;
   ncol = 0;
}

// constructor
Matrix::Matrix (const unsigned int nrow, const unsigned int ncol)
   :
   nrow (nrow),
   ncol (ncol)
{
   assert (nrow > 0);
   assert (ncol > 0);
   data = new double [nrow*ncol];
}

unsigned int Matrix::index (const unsigned int i, const unsigned int j)
{
   return i + nrow * j;
}

// assign one matrix to another
Matrix& Matrix::operator= (const Matrix& rhs)
{

   assert (nrow == rhs.nrow);
   assert (ncol == rhs.ncol);

   unsigned int n = nrow * ncol;
   for(unsigned int i=0; i<n; ++i)
      data[i] = rhs.data[i];
   
   return *this;
}

// assign one matrix to scalar
Matrix& Matrix::operator= (const double scalar)
{

   assert (nrow > 0);
   assert (ncol > 0);

   unsigned int n = nrow * ncol;
   for(unsigned int i=0; i<n; ++i)
      data[i] = scalar;
   
   return *this;
}

// add two matrices: result = this + mat2
Matrix Matrix::operator+ (const Matrix& mat2) const
{

   assert (nrow == mat2.nrow);
   assert (ncol == mat2.ncol);

   Matrix result (nrow, ncol);

   unsigned int n = nrow * ncol;
   for(unsigned int i=0; i<n; ++i)
      result.data[i] = data[i] + mat2.data[i];
   
   return result;
}

// add two matrices: this = this + mat2
Matrix& Matrix::operator+= (const Matrix& mat2)
{

   assert (nrow == mat2.nrow);
   assert (ncol == mat2.ncol);

   unsigned int n = nrow * ncol;
   for(unsigned int i=0; i<n; ++i)
      data[i] += mat2.data[i];
   
   return *this;
}

// subtract two matrices: result = this - mat2
Matrix Matrix::operator- (const Matrix& mat2) const
{

   assert (nrow == mat2.nrow);
   assert (ncol == mat2.ncol);

   Matrix result (nrow, ncol);

   unsigned int n = nrow * ncol;
   for(unsigned int i=0; i<n; ++i)
      result.data[i] = data[i] - mat2.data[i];
   
   return result;
}

// multiply two matrices element by element: result = this * mat2
Matrix Matrix::operator* (const Matrix& mat2) const
{

   assert (nrow == mat2.nrow);
   assert (ncol == mat2.ncol);

   Matrix result (nrow, ncol);

   unsigned int n = nrow * ncol;
   for(unsigned int i=0; i<n; ++i)
      result.data[i] = data[i] * mat2.data[i];
   
   return result;
}

// add two matrices: this = this - mat2
Matrix& Matrix::operator-= (const Matrix& mat2)
{

   assert (nrow == mat2.nrow);
   assert (ncol == mat2.ncol);

   unsigned int n = nrow * ncol;
   for(unsigned int i=0; i<n; ++i)
      data[i] -= mat2.data[i];
   
   return *this;
}

// multiply matrix with scalar: result = this * scalar
Matrix Matrix::operator* (const double scalar) const
{

   Matrix result (nrow, ncol);

   unsigned int n = nrow * ncol;
   for(unsigned int i=0; i<n; ++i)
      result.data[i] = scalar * data[i];
   
   return result;
}

// multiply matrix with scalar: this = this * scalar
Matrix& Matrix::operator*= (double scalar)
{

   unsigned int n = nrow * ncol;
   for(unsigned int i=0; i<n; ++i)
      data[i] *= scalar;
   
   return *this;
}

// destructor
Matrix::~Matrix ()
{
   if(nrow*ncol > 0)
      delete [] data;
}

// access matrix element (i,j)
double& Matrix::operator() (const unsigned int i, const unsigned int j) const
{
   return data[i + nrow*j];
}

// access matrix element (i,j)
double& Matrix::operator() (const unsigned int i, const unsigned int j)
{
   return data[i + nrow*j];
}

// allocate memory for matrix which has already been declared
void Matrix::allocate (const unsigned int ni, const unsigned int nj)
{
   assert (nrow == 0 && ncol == 0);
   assert (ni > 0 && nj > 0);

   nrow = ni;
   ncol = nj;
   data = new double[nrow*ncol];
}

// dot product of two matrices, element-by-element
double Matrix::dot (const Matrix &mat)
{
   assert (nrow == mat.nrow);
   assert (ncol == mat.ncol);

   double result = 0.0;
   unsigned int n = nrow * ncol;
   for(unsigned int i=0; i<n; ++i)
      result += data[i] * mat.data[i];

   return result;
}
