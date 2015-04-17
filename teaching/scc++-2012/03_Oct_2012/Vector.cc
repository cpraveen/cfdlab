#include <iostream>
#include <cassert>

#include "Vector.h"

//-----------------------------------------------------------------------------
// Constructor
//-----------------------------------------------------------------------------
template <class T>
Vector<T>::Vector (unsigned int nrow)
   :
      nrow (nrow)
{
   assert (nrow > 0);
   val.resize (nrow);
}

//-----------------------------------------------------------------------------
// Return size of vector
//-----------------------------------------------------------------------------
template <class T>
unsigned int Vector<T>::size () const
{
   return nrow;
}

//-----------------------------------------------------------------------------
// Return i'th value
//-----------------------------------------------------------------------------
template <class T>
T Vector<T>::operator() (unsigned int i) const
{
   return val[i];
}

//-----------------------------------------------------------------------------
// Return reference to i'th value
//-----------------------------------------------------------------------------
template <class T>
T& Vector<T>::operator() (unsigned int i)
{
   return val[i];
}

//-----------------------------------------------------------------------------
// Set vector to a constant value
//-----------------------------------------------------------------------------
template <class T>
Vector<T>& Vector<T>::operator= (const T scalar)
{
   for(unsigned int i=0; i<nrow; ++i)
      val[i] = scalar;
   return *this;
}

//-----------------------------------------------------------------------------
// Instantiation
//-----------------------------------------------------------------------------
template class Vector<int>;
template class Vector<double>;
