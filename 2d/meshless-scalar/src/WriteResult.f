c======================================================================
      subroutine WriteResult
c======================================================================
      implicit none
      include 'vars.h'

      integer i, vig

      vig = 10

      open(unit=vig, file="result.vig")
      write(vig,*)'points',np0
      do i=1,np0
         write(vig,*) coord(1,i), coord(2,i)
      enddo

      write(vig,*)'triangles',ntri
      do i=1,ntri
         write(vig,*) tri(1,i), tri(2,i), tri(3,i)
      enddo

      write(vig,*)'scalars  Sol'
      do i=1,np0
         write(vig,*) var(i)
      enddo
      write(vig,*)'end_block'

      close(vig)

      return
      end

c     write results for gnuplot
      subroutine WriteGNU
      implicit none
      include 'vars.h'

      integer gnu, i, j, c

      gnu = 30
      open(unit=gnu, file='result.dat')
      c = 0
      do i=1,nx
         do j=1,ny
            c = c + 1
            write(gnu,*) coord(1,c), coord(2,c), var(c)
         enddo
         write(gnu,*)
      enddo
      close(gnu)

      return
      end
