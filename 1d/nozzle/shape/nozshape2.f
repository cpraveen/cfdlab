C Set nozzle shape
      program main
      implicit none
      integer          nf
      parameter(nf=101)
      double precision x(nf), a(nf)

      integer          i
      double precision t, L, dx, ain, aout, da

      L   = 10.0d0
      ain = 0.21
      aout= 0.349

      print*,'Enter perturbation'
      read*,da

      dx  = L/(nf-1)
      open(unit=10, file='shape.dat')
      do i=1,nf
         x(i) = dx*(i-1)
         t    = x(i)/L
         a(i) = (-0.278d0+da)*t**3 + (0.417d0-da)*t**2 + 0.21d0
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
