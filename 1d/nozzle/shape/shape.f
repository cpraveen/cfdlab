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

      print*,'Enter control parameters c, d'
      read*,  dc, dd

      call shape(nf, L, ain, aout, dc, dd, a)

      stop
      end

      subroutine shape(nf, L, ain, aout, dc, dd, a)
      implicit none
      integer          nf
      double precision L, ain, aout, dc, dd, a(nf)

      integer          i
      double precision da, db, x(nf), dx

      db  = (aout - ain)/(dtanh(10.d0*dc - dd) - dtanh(-dd))
      da  = ain - db*dtanh(-dd)
      dx  = L/(nf-1)
      do i=1,nf
         x(i) = dx*(i-1)
         a(i) = da + db*dtanh(dc*x(i) - dd)
      enddo

      return
      end
