      subroutine check(xf,yf)
      use var
      implicit none
      real,intent(in) :: xf(1:imax,1:jmax), yf(1:imax,1:jmax)

      open(10,file="naca.plt")
      write(10,*)'TITLE = "NACA0012 2D structured grid"'
      write(10,*)'VARIABLES = "x", "y"'
      write(10,*)'ZONE STRANDID=1, SOLUTIONTIME=,t,',' I=',imax,&
                 ', J=',jmax,', DATAPACKING=POINT'
      do j=1,jmax
      do i=1,imax
         write(10,'(2E24.14)') xf(i,j), yf(i,j)
      enddo
      enddo
      close(10)
      end subroutine check
