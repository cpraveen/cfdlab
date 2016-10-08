      program converter
      use var
      implicit none
      
      call read_input

      write(*,*)'Reading mesh from file ',trim(ifilename)

      open(10,file=ifilename,status='old')
      read(10,*) nblocks
      read(10,*) imax, jmax
      allocate(x(1:imax, 1:jmax), y(1:imax, 1:jmax))
      do k=1,nblocks
         read(10,*) ((x(i,j),i=1,imax),j=1,jmax), &
                    ((y(i,j),i=1,imax),j=1,jmax)
      enddo
      close(10)

      write(*,'(X,A12,2I6)') 'imax, jmax =', imax, jmax

      call check(x,y)
      call msh_ascii(x,y)
      end program
