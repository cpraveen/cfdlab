      subroutine msh_ascii(xg,yg)
      use var
      implicit none
      integer :: kk, nnodes, nelems, eltype
      real, intent(in) :: xg(1:imax,1:jmax), yg(1:imax,1:jmax)
      integer, dimension(:,:), allocatable :: ind
      integer, dimension(:), allocatable :: tag

      write(*,'(X,A12,I5)') 'TE start at ',wst
      write(*,'(X,A12,I5)') 'TE end   at ',wen
      write(*,*)'Writing mesh into ',trim(ofilename)

      nnodes = imax*jmax - wst
      nelems = (imax-1)*(jmax-1) + (wen-wst) + (imax-1)   +  2*(jmax-1)
      !-------quadrangles(1)-------wall(2)----farfield(3)----farfield(4)-----

      allocate(tag(1:ntags),ind(imax+1,jmax+1))

      open(8, file=ofilename)
      write(8,'(A11)') "$MeshFormat"
      write(8,'(F4.2,X,I1,X,I1)') version, file_type, data_size
      write(8,'(A14)') "$EndMeshFormat"

      write(8,'(A6)') "$Nodes"
      write(8,'(I6)') nnodes
      k = 1
      j = 1
       do i = 1,wen-1
          write(8,'(I6,X,E20.14,2X,E20.14,3X,F3.1)') k,xg(i,j),yg(i,j),z
          ind(i,j) = k
          k = k+1
       enddo
       do j = 2,jmax
       do i = 1,imax
          write(8,'(I6,X,E22.16,2X,E22.16,3X,F3.1)') k,xg(i,j),yg(i,j),z
          ind(i,j) = k
          k = k+1
       enddo
       enddo
      write(8,'(A9)') "$EndNodes"

      write(8,'(A9)') "$Elements"
      write(8,'(I6)') nelems

      eltype = 3 ; tag(1) = 1 ; tag(2) = 1  ! for grid
      k = 1
      j = 1    
      do i = 1,imax-1   
         if(i .lt. wen-1) then
           write(8,'(9(I6,X))') k,eltype,ntags,tag(1),tag(2),ind(i,j), &
                                ind(i+1,j),ind(i+1,j+1),ind(i,j+1)
           k = k+1
         else
           write(8,'(9(I6,X))') k,eltype,ntags,tag(1),tag(2), &
                                ind(wen-i+wst,j),ind(wen-i+wst-1,j), &
                                ind(i+1,j+1),ind(i,j+1)
           k = k+1
         endif
      enddo
      do j = 2,jmax-1
      do i = 1,imax-1
         write(8,'(9(I6,X))') k,eltype,ntags,tag(1),tag(2),ind(i,j), &
                              ind(i+1,j),ind(i+1,j+1),ind(i,j+1)
         k = k+1
      enddo
      enddo

      eltype = 1 ; tag(1) = 2 ; tag(2) = 2  ! for airfoil wall
      j = 1
      do i = wst, wen-2
       write(8,'(7(I6,X))') k,eltype,ntags,tag(1),tag(2),ind(i,j),&
                            ind(i+1,j)
       k = k+1
      enddo
      i = wen-1
      write(8,'(7(I6,X))') k,eltype,ntags,tag(1),tag(2),ind(i,j),ind(wst,j)
      k = k+1

      eltype = 1 ; tag(1) = 3 ; tag(2) = 3  ! for farfield_1 (curved)
      j = jmax
      do i = 1, imax-1
       write(8,'(7(I6,X))') k,eltype,ntags,tag(1),tag(2),ind(i,j), &
                            ind(i+1,j)
       k = k+1
      enddo

      eltype = 1 ; tag(1) = 4 ; tag(2) = 4  ! for farfield_2 (vertical at bottom)
      i = 1
      do j = 1, jmax-1
       write(8,'(7(I6,X))') k,eltype,ntags,tag(1),tag(2),ind(i,j), &
                            ind(i,j+1)
       k = k+1
      enddo

      eltype = 1 ; tag(1) = 4 ; tag(2) = 4  ! for farfield_3 (vertical at upside)
      i = imax
      do j = 1, jmax-1
       if(j .eq. 1)then
          write(8,'(7(I6,X))') k,eltype,ntags,tag(1),tag(2),ind(1,j), &
                               ind(i,j+1)
       else
          write(8,'(7(I6,X))') k,eltype,ntags,tag(1),tag(2),ind(i,j), &
                               ind(i,j+1)
       endif
       k = k+1
      enddo

      write(8,'(A12)') "$EndElements"
      close(8)

      end subroutine msh_ascii
