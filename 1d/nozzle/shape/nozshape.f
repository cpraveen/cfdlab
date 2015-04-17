C Set nozzle shape
      program main
      implicit none
      integer          nf
      parameter(nf=101)
      double precision x(nf), a(nf)

      integer          i
      double precision L, dx, ain, aout, da, db, dc, dd

      L   = 10.0d0
      ain = 1.0512d0
      aout= 1.75d0

c     Target shape
      dc  = 0.8d0
      dd  = 4.0d0
c     Start shape: try dc=1.0, dd=3.8

      print*,'Enter control parameters c, d'
      read*,  dc, dd

      db  = (aout - ain)/(dtanh(L*dc - dd) - dtanh(-dd))
      da  = ain - db*dtanh(-dd)
      dx  = L/(nf-1)
      open(unit=10, file='shape.dat')
      do i=1,nf
         x(i) = dx*(i-1)
         a(i) = da + db*dtanh(dc*x(i) - dd)
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
      print*,'b                  =',db
      print*,'c                  =',dc
      print*,'d                  =',dd

      print*,'Shape written into shape.dat'

      stop
      end
