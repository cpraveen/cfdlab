      subroutine InterpolateFlux(x1, x2, z1, z2, xgauss, flux)
      implicit none
      include 'param.inc'
      double precision x1, x2, z1(*), z2(*), xgauss, flux(*)

      integer i
      double precision z(3), rho, u, p, e, hnozzle
      double precision nozzle

      do i=1,3
         z(i) = z1(i) + (xgauss - x1)*(z2(i)-z1(i))/(x2-x1)
      enddo

      hnozzle = nozzle(xgauss)

      rho = z(1)**2
      u = z(2)/z(1)
      p = gam1*(z(1)*z(3) - 0.5d0*z(2)**2)/gam
      e = p/gam1 + 0.5d0*rho*u**2
      flux(1) = hnozzle*z(1)*z(2)
      flux(2) = hnozzle*(p + z(2)**2)
      flux(3) = hnozzle*(e + p)*u

      end
