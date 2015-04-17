#include <iostream>
#include <cassert>
#include <vector>

#include "sparse_matrix.h"
#include "Vector.h"

//-----------------------------------------------------------------------------
// Constructor
//-----------------------------------------------------------------------------
template <class T>
SparseMatrix<T>::SparseMatrix (std::vector<unsigned int>& row_ptr, 
                               std::vector<unsigned int>& col_ind, 
                                          std::vector<T>& val)
   :
   nrow (row_ptr.size()-1),
   row_ptr (row_ptr),
   col_ind (col_ind),
   val (val)
{
   assert (row_ptr.size() >= 2);
   assert (col_ind.size() > 0);
   assert (col_ind.size() == val.size());
   assert (row_ptr[nrow] == val.size());

   for(unsigned int i=0; i<col_ind.size(); ++i)
      assert (col_ind[i] >= 0 && col_ind[i] < nrow);
}

//-----------------------------------------------------------------------------
// Get element value of A(i,j)
//-----------------------------------------------------------------------------
template <class T>
T SparseMatrix<T>::operator()(unsigned int i, 
                              unsigned int j) const
{
   unsigned int row_beg = row_ptr[i];
   unsigned int row_end = row_ptr[i+1];
   for(unsigned int d=row_beg; d<row_end; ++d)
      if(col_ind[d] == j) 
         return val[d];
   return 0;
}

//-----------------------------------------------------------------------------
// y = A*x
//-----------------------------------------------------------------------------
template <class T>
void SparseMatrix<T>::multiply(const Vector<T>& x, 
                                     Vector<T>& y) const
{
   assert (x.size() == nrow);
   assert (x.size() == y.size());
   for(unsigned int i=0; i<nrow; ++i)
   {
      y(i) = 0.0;
      unsigned int row_beg = row_ptr[i];
      unsigned int row_end = row_ptr[i+1];
      for(unsigned int j=row_beg; j<row_end; ++j)
         y(i) += val[j] * x(col_ind[j]);
   }
}

//-----------------------------------------------------------------------------
// Instantiation
//-----------------------------------------------------------------------------
template class SparseMatrix<int>;
template class SparseMatrix<double>;
