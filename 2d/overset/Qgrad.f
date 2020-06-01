c======================================================================
      subroutine Qgrad(n, x0, y0, Q0, x1, y1, Q1, Qx, Qy)
c======================================================================
c n                 = Number of neighbours
c x0,y0             = coordinates of point at which grad is required
c Q0(nvar)          = primitive variables
c x1(n),y1(n)       = coordinates of neighbouring points
c Q1(nvar,n)        = prim vars at neighbouring points
c Qx(nvar),Qy(nvar) = Grad at (x0,y0)
c======================================================================
      implicit none
      include 'misc.h'
      integer n
      real    x0, y0, Q0(*), x1(*), y1(*), Q1(nvar,*), Qx(*), Qy(*)

c     local variables
      real    cx(n), cy(n)

c     compute coefficients in least squares formula for derivative
      call LSCoeff(n, x0, y0, x1, y1, cx, cy)
      call Qgrad2(n, Q0, Q1, cx, cy, Qx, Qy)

      return
      end

c======================================================================
      subroutine Qgrad2(n, Q0, Q1, cx, cy, Qx, Qy)
c======================================================================
c Given states and coefficients, computes derivatives
c======================================================================
      implicit none
      include 'misc.h'
      integer n
      real    Q0(*), Q1(nvar,*), cx(*), cy(*), Qx(*), Qy(*)

      integer i, j
      real    dQ

c     compute derivative
      do i=1,nvar
         Qx(i) = 0.0
         Qy(i) = 0.0
         do j=1,n
            dQ    = Q1(i,j) - Q0(i)
            Qx(i) = Qx(i) + cx(j)*dQ
            Qy(i) = Qy(i) + cy(j)*dQ
         enddo
      enddo

      return
      end
