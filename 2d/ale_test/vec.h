#ifndef __VEC_H__
#define __VEC_H__

#include <cmath>

// 3-D vector
class Vector
{
   public:
      double x, y;
      Vector& operator=  (const Vector rhs);
      Vector& operator=  (const double scalar);
      Vector& operator+= (const Vector rhs);
      Vector& operator-= (const Vector rhs);
      Vector& operator*= (const double scalar);
      Vector& operator/= (const double scalar);
      Vector  operator/  (const double scalar) const;
      Vector  operator*  (const double scalar) const;
      double  operator*  (const Vector vec) const; // Dot product of two vectors
      Vector  operator+  (const Vector vec) const;
      Vector  operator-  (const Vector vec) const;
      double  operator^  (const Vector vec) const; // Cross product of two vectors
      double  square () const;
      double  norm () const;
};

// Assign one vector to another
inline
Vector& Vector::operator= (const Vector rhs){
   x = rhs.x;
   y = rhs.y;

   return *this;
}

// Assign one vector to another
inline
Vector& Vector::operator= (const double scalar){
   x = y = scalar;

   return *this;
}

// Add vector to given vector: this = this + rhs
inline
Vector& Vector::operator+= (const Vector rhs){
   x += rhs.x;
   y += rhs.y;

   return *this;
}

// Subtract vector from given vector: this = this - rhs
inline
Vector& Vector::operator-= (const Vector rhs){
   x -= rhs.x;
   y -= rhs.y;

   return *this;
}

// Multiply vector by scalar and copy result to same vector
inline
Vector& Vector::operator*= (const double scalar){
   x *= scalar;
   y *= scalar;

   return *this;
}

// Divide vector by scalar and copy result to same vector
inline
Vector& Vector::operator/= (const double scalar){
   x /= scalar;
   y /= scalar;

   return *this;
}

// Add two vectors
inline
Vector Vector::operator+  (const Vector vec) const
{
   Vector result;

   result.x = x + vec.x;
   result.y = y + vec.y;

   return result;
}

// Subtract two vectors
inline
Vector Vector::operator-  (const Vector vec) const
{
   Vector result;

   result.x = x - vec.x;
   result.y = y - vec.y;

   return result;
}

// Divide a vector by a scalar
inline
Vector Vector::operator/ (const double scalar) const
{
   Vector result;
   result.x = x / scalar;
   result.y = y / scalar;

   return result;
}

// Multiply a vector by a scalar
inline
Vector Vector::operator* (const double scalar) const
{
   Vector result;
   result.x = x * scalar;
   result.y = y * scalar;

   return result;
}

// L2 norm square of vector
inline
double Vector::square () const
{
   return x*x + y*y;
}

// L2 norm of vector
inline
double Vector::norm () const
{
   return std::sqrt(x*x + y*y);
}

// Dot product of two vectors
inline
double Vector::operator* (const Vector vec) const
{
   return x * vec.x + y * vec.y;
}

// Cross product of two vectors
inline
double Vector::operator^ (const Vector vec) const
{
   return  x * vec.y - y * vec.x;
}

#endif
