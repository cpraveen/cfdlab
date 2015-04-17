c     Generate overset grids in 1-d !!!
      program grid
      implicit none
      include 'param.h'

      integer ifid, i, j, ig1, ig2
      integer irp, irg, idg, id1, id2
      real    ds, fact


c     call SetGrid1
c     call SetGrid2
c     call SetGrid3
      call SetGrid4

c     Generate grid
      do i=1,ngrid
         if(np(i).gt.npmax)then
            print*,'Increase npmax to atleast ',np(i)
            stop
         endif
         dx(i) = (xmax(i) - xmin(i))/(np(i) - 1)
         do j=1,np(i)
            x(i,j) = xmin(i) + (j-1)*dx(i)
         enddo
         print*,'grid=',i,' dx=',dx(i)
      enddo

c     Set update type
      utype(:,:)        = UPDATE
      utype(1, 1      ) = FIXED
      utype(1, 2      ) = FIXED
      utype(1, np(1)  ) = OVRLAP
      utype(1, np(1)-1) = OVRLAP
      utype(2, 1      ) = OVRLAP
      utype(2, 2      ) = OVRLAP
      utype(2, np(2)  ) = FIXED
      utype(2, np(2)-1) = FIXED

c     Write out the grid
      ifid    = 10
      open(ifid, file='grids.dat')
      write(ifid,*) ngrid
      write(ifid,*) (np(i), i=1,ngrid)
      write(ifid,*) (dx(i), i=1,ngrid)
      do j=1,ngrid
         do i=1,np(j)
            write(ifid,*) x(j,i), utype(j,i)
         enddo
      enddo
      close(ifid)

      open(ifid, file='g.dat')
      do j=1,ngrid
         do i=1,np(j)
            write(ifid,*) x(j,i)
         enddo
         write(ifid,*)
      enddo
      close(ifid)

c     Count number of overlap points
      nopts = 0

      if(ngrid.eq.1) goto 100

c     Find overlap information
      ig1 = 1
      ig2 = 2
      do j=1,np(ig1)
         if(utype(ig1,j).eq.OVRLAP)then
            nopts = nopts + 1
            do i=1,np(ig2)
               if(x(ig2,i).le.x(ig1,j) .and. x(ig1,j).lt.x(ig2,i+1))then
                  olap(1,nopts) = ig1
                  olap(2,nopts) = j
                  olap(3,nopts) = ig2
                  olap(4,nopts) = i
                  olap(5,nopts) = i+1
                  olap(6,nopts) = i-1
                  olap(7,nopts) = i+2
                  print*,ig1,j,ig2,i,i+1
                  print*, x(ig1,j), x(ig2,i), x(ig2,i+1)
               endif
            enddo
         endif
      enddo

      ig1 = 2
      ig2 = 1
      do j=1,np(ig1)
         if(utype(ig1,j).eq.OVRLAP)then
            nopts = nopts + 1
            do i=1,np(ig2)
               if(x(ig2,i).le.x(ig1,j) .and. x(ig1,j).lt.x(ig2,i+1))then
                  olap(1,nopts) = ig1
                  olap(2,nopts) = j
                  olap(3,nopts) = ig2
                  olap(4,nopts) = i
                  olap(5,nopts) = i+1
                  olap(6,nopts) = i-1
                  olap(7,nopts) = i+2
                  print*,ig1,j,ig2,i,i+1
                  print*, x(ig1,j), x(ig2,i), x(ig2,i+1)
               endif
            enddo
         endif
      enddo
      print*,'Number of overlap points =',nopts


c     Find overlap factors
      do i=1,nopts
         irg = olap(1,i)
         irp = olap(2,i)

         idg = olap(3,i)
         id1 = olap(4,i)
         id2 = olap(5,i)
         ds  =   x(idg,id2) - x(idg,id1)
         fact= ( x(irg,irp) - x(idg,id1) )/ds
         print*,i,fact
      enddo

100   continue

      open(20, file='overlap.dat')
      write(20,*) nopts
      do i=1,nopts
         write(20,'(7I6)') (olap(j,i), j=1,7)
      enddo
      close(20)


      stop
      end
c------------------------------------------------------------------------------
c------------------------------------------------------------------------------
      subroutine SetGrid1
      implicit none
      include 'param.h'

c     Number of grids
      ngrid   = 1

c     Grid 1
      xmin(1) = 0.0
      xmax(1) = 2.0
      np  (1) = 201

      np  (2) = 2

      return
      end
c------------------------------------------------------------------------------
c------------------------------------------------------------------------------
      subroutine SetGrid2
      implicit none
      include 'param.h'

c     Number of grids
      ngrid   = 2

c     Grid 1
      xmin(1) = 0.0
      xmax(1) = 1.0
      np  (1) = 51

c     Grid 2
      xmin(2) = 0.75
      xmax(2) = 2.01
      np  (2) = 64

      return
      end
c------------------------------------------------------------------------------
c------------------------------------------------------------------------------
      subroutine SetGrid3
      implicit none
      include 'param.h'

c     Number of grids
      ngrid   = 2

c     Grid 1
      xmin(1) = 0.0
      xmax(1) = 1.0
      np  (1) = 51

c     Grid 2
      xmin(2) = 0.76
      xmax(2) = 2.0
      np  (2) = 63

      return
      end
c------------------------------------------------------------------------------
      subroutine SetGrid4
      implicit none
      include 'param.h'

c     Number of grids
      ngrid   = 2

c     Grid 1
      xmin(1) = 0.0
      xmax(1) = 1.0
      np  (1) = 101

c     Grid 2
      xmin(2) = 0.755
      xmax(2) = 2.0
      np  (2) = 126

      return
      end
