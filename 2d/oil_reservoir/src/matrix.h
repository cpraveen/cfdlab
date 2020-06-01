#ifndef __MATRIX_H__
#define __MATRIX_H__

class Matrix
{
   public:
      Matrix ();
      Matrix (const unsigned int nrow, const unsigned ncol);
      Matrix& operator= (const Matrix&);
      Matrix& operator= (const double);
      Matrix  operator+ (const Matrix&) const; // add two matrices
      Matrix  operator- (const Matrix&) const; // subtract two matrices
      Matrix  operator* (const Matrix&) const; // element-by-element multiply
      Matrix& operator*= (double);             // multiply by scalar
      Matrix& operator+= (const Matrix&);      // add a matrix
      Matrix& operator-= (const Matrix&);      // subtract a matrix
      Matrix  operator* (double) const;        // multiply by scalar
      double dot (const Matrix&);
      ~Matrix ();
      double& operator() (const unsigned int i, const unsigned int j) const;
      double& operator() (const unsigned int i, const unsigned int j);
      void allocate (const unsigned int, const unsigned int);
      unsigned int index (const unsigned int, const unsigned int);

   private:
      unsigned int nrow, ncol;
      double* data;
};

#endif
