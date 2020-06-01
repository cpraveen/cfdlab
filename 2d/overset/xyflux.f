c======================================================================
      subroutine xyflux(q, fx, fy)
c======================================================================
c Computes fluxes in x and y direction given
c q = (rho, u, v, p)
c Requires gm1 = gamma - 1.0, defined locally. But replace it later.
c======================================================================
      implicit none
      real q(*), fx(*), fy(*)

      real E, gm1

      gm1 = 1.4 - 1.0

c     Energy per unit volume
      E     = q(4)/gm1 + 0.5*q(1)*(q(2)**2 + q(3)**2)

c     x flux
      fx(1) = q(1)*q(2)
      fx(2) = q(4) + q(1)*q(2)**2
      fx(3) = q(1)*q(2)*q(3)
      fx(4) = (q(4) + E)*q(2)

c     y flux
      fy(1) = q(1)*q(3)
      fy(2) = q(1)*q(2)*q(3)
      fy(3) = q(4) + q(1)*q(3)**2
      fy(4) = (q(4) + E)*q(3)

      return
      end
