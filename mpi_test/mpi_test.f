C
C     Driver routine to run flow solver for different ALPHA (AOA)
C
C     Specify min, max and increment for the AOA
C     Cl and Cd will be written into alpha.dat
C
      program alpha_sweep
      implicit none
      include 'mpif.h'

      integer :: nproc, rank, ierr, i
      real    :: x(1000)

      call MPI_init(ierr)
      call MPI_comm_size(MPI_comm_world, nproc, ierr)
      call MPI_comm_rank(MPI_comm_world, rank, ierr)

      call system('hostname')

      do i=1,nproc
         x(i) = 0.0
      enddo
      x(rank+1) = rank

      do i = 1, nproc
          call MPI_Bcast(x(i), 1, MPI_real, rank, MPI_comm_world
     +      ,ierr)
      end do
      
      
      call MPI_finalize(ierr)

      stop
      end program
