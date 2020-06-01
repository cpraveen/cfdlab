c======================================================================
c Set initial condition
      subroutine SetInitCond
c======================================================================
      implicit none
      include 'vars.h'

      integer          i, p
      double precision x, y, InitCond

      print*,'Setting initial condition ...'

      varmin0 = 1.0d20
      varmax0 =-1.0d20
      do i=1,np0
         x       = coord(1,i)
         y       = coord(2,i)
         var(i)  = InitCond(x, y)
         varmin0 = dmin1(varmin0, var(i))
         varmax0 = dmax1(varmax0, var(i))
      enddo

c     periodic condition
      do i=np0+1,npts
         p = conn(1,i)
         var(i) = var(p)
      enddo

      print*,'Initial min,max =',varmin0,varmax0

      return
      end

c======================================================================
c Initial condition
      double precision function InitCond(x, y)
c======================================================================
      implicit none
      include 'vars.h'
      double precision x, y

      double precision r, r2

      if(ictype.eq.1)then
c        use xmin=0, xmax=1, ymin=0, ymax=1
         r2 = (x-0.5d0)**2 + (y-0.5d0)**2
         r  = dsqrt(r2)
         if(r.le.0.25d0)then
            InitCond = (0.25d0**2 - r2)/0.25d0**2
            if(InitCond.lt.0.0d0)InitCond = 0.0d0
         else
            InitCond = 0.0d0
         endif
      else
         print*,'ictype is not known, ictype =',ictype
         stop
      endif

      return
      end
