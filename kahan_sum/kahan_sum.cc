// Kahan summation to avoid roundoff errors, see
// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
// For a fortran example, see
// http://ossanworld.com/cfdbooks/cfdcodes/kahan_sum_example/kahan_sum_example.f90
#include <iostream>
#include <iomanip>
#include <vector>
#include "kahan_sum.h"

using namespace std;

// Compute 1 + sum(i=0,n-1) v[i]
template <typename T>
void test()
{
   unsigned int n = 100000000;
   T one = 1.0;
   T a = 1.234567891234567e-2;
   vector<T> v(n, a);

   KahanSum<T> ksum(one); // Kahan sum
   T sum = one;           // Ordinary sum
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
