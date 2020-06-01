c     find initial global time step allowed by stability condition
      subroutine InitTimeStep
      implicit none
      include 'vars.h'

      integer          i, j
      double precision x, y, vx, vy, v, vsum, dtlocal, dtglobal

      dtglobal = 1.0d20

      do i=1,np0
         vsum = 0.0d0
         x    = coord(1,i)
         y    = coord(2,i)
         call AdvecVel(x,y,vx,vy)
         do j=1,nnbr(i)
            v    = ax(j,i)*vx + ay(j,i)
            vsum = vsum + dabs(v)
         enddo
         vsum    = vsum + 1.0d-16
         dtlocal = 1.0d0/vsum
         dtglobal= dmin1(dtglobal, dtlocal)
      enddo

      print*,'Chosen time step =',dt
      print*,'CFL time step    =',dtglobal

      if(dt.gt.dtglobal+1.0d-15)then
         stop 'Time step is greater than stability limit'
      endif

      return
      end
