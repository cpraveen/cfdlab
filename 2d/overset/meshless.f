c======================================================================
      subroutine meshless(dt, n, x0, y0, Q0, x1, y1, Q1, Qx1, Qy1)
c======================================================================
c n = number of points req meshless update
c dt(n) = local time-step
c x0(n),y0(n) = coordinates of pts req meshless update
c x1(nnbr,n),y1(nnbr,n) = coordinates of neighbouring points
c Q0(nvar,n)
c Q1(nvar,nnbr,n)
c Qx1, Qy1(nvar,nnbr,n) = derivatives at neighbouring pts
c======================================================================
      implicit none
      include 'misc.h'
      integer n
      real    dt(*), x0(*), y0(*), Q0(nvar,*), x1(nnbr,*), y1(nnbr,*),
     +        Q1(nvar,nnbr,*), Qx1(nvar,nnbr,*), Qy1(nvar,nnbr,*)

c     local variables
      integer i, j, k
      real    cx(25), cy(25), Qx0(nvar), Qy0(nvar), Ql(nvar), Qr(nvar),
     +        fx0(nvar), fy0(nvar), flr(nvar), residu(nvar)

c     for each meshless point
      do i=1,n

c        find coefficients
         call LSCoeff(nnbr, x0(i), y0(i), x1(1,i), y1(1,i), cx, cy)

c        find q-derivatives
         call Qgrad2(nnbr, Q0(1,i), Q1(1,1,i), cx, cy, Qx0, Qy0)

c        find xyflux
         call xyflux(Q0(1,i), fx0, fy0)

         do k=1,nvar
            residu(k) = 0.0
         enddo

c        for each neighbour
         do j=1,nnbr

c           find left-right states
            call reconstruct(x0(i), y0(i), Q0(1,i), Qx0, Qy0, 
     +                       x1(j,i), y1(j,i), Q1(1,j,i), Qx1(1,j,i),
     +                       Qy1(1,j,i), Ql, Qr)

c           find flux
            call roeflx(flr, Ql, Qr, cx(j), cy(j), 0.0)

c           update residue
            do k=1,nvar
               residu(k) = residu(k) + 2.0*(flr(k) - cx(j)*fx0(k) -
     +                                               cy(j)*fy0(k))
            enddo

         enddo

c        update solution
         do k=1,nvar
            Q0(k,i) = Q0(k,i) - dt(i)*residu(k)
         enddo

      enddo

      return
      end
