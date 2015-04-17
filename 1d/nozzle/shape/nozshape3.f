C Set nozzle shape
      program main
      implicit none
      integer          nf
      parameter(nf=101)
      double precision x(nf), a(nf)

      integer          i
      double precision t, L, dx, ain, aout, da

      L   = 1.0d0
      ain = 1.0304d0
      aout= 1.1022d0

      print*,'Enter perturbation'
      read*,da

      dx  = L/(nf-1)
      open(unit=10, file='shape.dat')
      do i=1,nf
         x(i) = dx*(i-1)
         t    = x(i)
         a(i) = ain + 0.0718*t**3*(10.0d0-15.0d0*t+6.0d0*t**2)
     +          - 0.9482d0*t**3*(1.0d0-t)**3
         write(10,'(2e20.10)') x(i), a(i)
      enddo
      aout= a(nf)
      close(10)

      print*,'Number of points   =',nf
      print*,'Length of nozzle   =',L
      print*,'dx                 =',dx
      print*,'Inlet area         =',ain
      print*,'Outlet area        =',aout
      print*,'Area ratio         =',aout/ain

      print*,'Shape written into shape.dat'

      stop
      end
