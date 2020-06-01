      subroutine Gradients
      implicit none
      include 'vars.h'

      integer i, j, p

      do i=1,npts
         if(ptype(i).eq.INTERIOR)then
            varx(i) = 0.0d0
            vary(i) = 0.0d0
            do j=1,nnbr(i)
               p       = conn(j,i)
               varx(i) = varx(i) + ax(j,i)*(var(p) - var(i))
               vary(i) = vary(i) + ay(j,i)*(var(p) - var(i))
            enddo
         endif
      enddo

c     periodic condition
      do i=1,npts
         if(ptype(i).eq.PERIODIC)then
            p = conn(1,i)
            varx(i) = varx(p)
            vary(i) = vary(p)
         endif
      enddo

      return
      end
