// Kahan summation to avoid roundoff errors, see
// https://en.wikipedia.org/wiki/Kahan_summation_algorithm

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
