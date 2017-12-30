// Kahan summation to avoid roundoff errors, see
// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
// For a fortran example, see
// http://ossanworld.com/cfdbooks//cfdcodes/kahan_sum_example/kahan_sum_example.f90
#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

template <typename T>
class KahanSum
{
   public:
      KahanSum (T init=0.0);
      KahanSum& operator=  (const T &v);
      KahanSum& operator+= (const T &v);
      T result() const;

   private:
      T c;
      T sum;
};

template <typename T>
KahanSum<T>::KahanSum(T init)
{
   c   = 0.0;
   sum = init;
}

template <typename T>
KahanSum<T>& KahanSum<T>::operator+=(const T &v)
{
   T temp = v - c;
   T t    = sum + temp;
   c      = (t - sum) - temp;
   sum    = t;
   return *this;
}

template <typename T>
KahanSum<T>& KahanSum<T>::operator=(const T &v)
{
   c      = 0.0;
   sum    = v;
   return *this;
}

template <typename T>
T KahanSum<T>::result() const
{
   return sum;
}

// Compute 1 + sum(i=0,n-1) v[i]
template <typename T>
void test()
{
   unsigned int n = 100000000;
   T one = 1.0;
   T a = 1.234567891234567e-2;
   vector<T> v(n, a);

   KahanSum<T> ksum(one);
   T sum = one;
   for(size_t i=0; i<n; ++i)
   {
      ksum += v[i];
      sum  += v[i];
   }

   cout << setprecision(14);
   cout << "      sum = " << sum << endl;
   cout << "Kahan sum = " << ksum.result() << endl;
   cout << "Exact     = " << one + a*n << endl; 
}

int main()
{
   cout << "Single precision:" << endl;
   test<float>();

   cout << "Double precision:" << endl;
   test<double>();
}
