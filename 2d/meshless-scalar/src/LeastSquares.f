      subroutine LeastSquares
      implicit none
      include 'vars.h'

      integer          i, j, nn
      double precision x0, y0, x1(nnbrmax), y1(nnbrmax)

      do i=1,npts
         if(ptype(i).ne.INTERIOR) goto 100
         x0 = coord(1,i)
         y0 = coord(2,i)
         do j=1,nnbr(i)
            nn    = conn(j,i)
            x1(j) = coord(1,nn)
            y1(j) = coord(2,nn)
         enddo
         call LSCoeff(nnbr(i), x0, y0, x1, y1, ax(1,i), ay(1,i))
100      continue
      enddo

      return
      end
