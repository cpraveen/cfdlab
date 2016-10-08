      module var
      implicit none
      integer :: i, j, k, nblocks
      integer :: imax, jmax, wst, wen
      integer :: nphys, ntags 
      integer :: file_type, data_size, nperiodic

      real :: version
      real,parameter :: z = 0.0
      real, dimension(:,:), allocatable :: x,y

      character (len=40) :: ifilename, ofilename
      end module var
