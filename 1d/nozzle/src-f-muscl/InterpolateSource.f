      subroutine InterpolateSource(x1, x2, z1, z2, xgauss, advar,source)
      implicit none
      include 'param.inc'
      double precision x1, x2, z1(*), z2(*), xgauss, advar(*), source

      integer i
      double precision z(3), p, xgaussb

      do i=1,3
         z(i) = z1(i) + (xgauss - x1)*(z2(i)-z1(i))/(x2-x1)
      enddo

      p = gam1*(z(1)*z(3) - 0.5d0*z(2)**2)/gam

      call nozzle_x(xgauss, xgaussb, 1.0d0)

      source = advar(2) * xgaussb * p

      end
