c======================================================================
      real function valimit(a, b)
c======================================================================
c van-albada limiter function
c======================================================================
      implicit none
      real    a, b

      real    num, den, epsi
      parameter(epsi=1.0e-10)

      num     = 2.0*a*b + epsi
      den     = a**2 + b**2 + epsi
      valimit = max(0.0, num/den)

      return
      end
