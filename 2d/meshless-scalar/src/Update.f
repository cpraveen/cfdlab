      subroutine Update(rks)
      implicit none
      include 'vars.h'
      integer          rks

      integer          i
      double precision aa(10), bb(10)

c     for 3-stage scheme
      aa(1) = 0.0d0
      aa(2) = 3.0d0/4.0d0
      aa(3) = 1.0d0/3.0d0
      bb(1) = 1.0d0
      bb(2) = 1.0d0/4.0d0
      bb(3) = 2.0d0/3.0d0

      do i=1,npts
         if(ptype(i).eq.INTERIOR)then
            var(i) = aa(rks)*var_old(i) + bb(rks)*(var(i) - dt*res(i))
         endif
      enddo

c     apply periodic condition
      do i=1,npts
         if(ptype(i).eq.PERIODIC)then
            var(i) = var(conn(1,i))
         endif
      enddo

      return
      end
