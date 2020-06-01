c======================================================================
      subroutine reconstruct(x0, y0, Q0, Qx0, Qy0, x1, y1, Q1, Qx1, Qy1,
     +                       Ql, Qr)
c======================================================================
c MUSCL-type linear reconstruction with van-albada limiter
c Check for the constant kkk in the reconstruction. 
c======================================================================
      implicit none
      include 'misc.h'
      real    x0, y0, Q0(*), Qx0(*), Qy0(*), x1, y1, Q1(*), 
     +        Qx1(*), Qy1(*), Ql(*), Qr(*)

      integer i
      real    dx, dy, dQ, dm, dp, sl, sr, VALIMIT, kkk
      external VALIMIT

      kkk = 1.0/3.0

      dx    = x1 - x0
      dy    = y1 - y0

      do i=1,nvar
         dQ    = Q1(i) - Q0(i)
         dm    = 2.0*(dx*Qx0(i) + dy*Qy0(i)) - dQ
         dp    = 2.0*(dx*Qx1(i) + dy*Qy1(i)) - dQ
         sl    = VALIMIT(dm, dQ)
         sr    = VALIMIT(dp, dQ)
         Ql(i) = Q0(i) + 0.25*sl*( (1.0-kkk*sl)*dm + (1.0+kkk*sl)*dQ )
         Qr(i) = Q1(i) - 0.25*sr*( (1.0-kkk*sr)*dp + (1.0+kkk*sr)*dQ )
      enddo

      return
      end
