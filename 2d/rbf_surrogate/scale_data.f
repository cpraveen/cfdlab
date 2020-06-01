C-----------------------------------------------------------------------------
C Scale coordinates and function values
C NOTE: If scalex=0, then xmin and xmax are not computed.
C-----------------------------------------------------------------------------
      subroutine scale_xdata(scalex, nvar, npt, xs, xmin, xmax)
      implicit none
      include 'Param.h'
      integer          nvar, npt, scalex
      double precision xs(nvar,npt), xmin(nvar), xmax(nvar)

      integer          ivar, ipt
      double precision xrange
      
c Find bounding box in design space and scale
      do ivar=1,nvar
         xmin(ivar) = +1.0d20
         xmax(ivar) = -1.0d20
         do ipt=1,npt
            if(xs(ivar,ipt) .lt. xmin(ivar)) xmin(ivar) = xs(ivar,ipt)
            if(xs(ivar,ipt) .gt. xmax(ivar)) xmax(ivar) = xs(ivar,ipt)
         enddo
      enddo

      do ivar=1,nvar
         xrange = xmax(ivar) - xmin(ivar)
         if(xrange.lt.MACHEPS)then
            write(*,*)'###### Scale_data ######'
            write(*,*)' No variation in variable=',ivar
         endif
      enddo

c If coordinate values dont need to be scaled
      if(scalex .eq. 0)then
         print*,'Coordinates will not be scaled'
         return
      endif

c Scale the coordinates
      print*,'Scaling coordinates ...'
      do ivar=1,nvar
         xrange = xmax(ivar) - xmin(ivar)
         if(xrange .lt. MACHEPS) then
            write(*,*) '******* WARNING *******'
            write(*,*) 'Variable number ',ivar,' is constant'
            xrange = 1.0d0
         endif
         xrange = 1.0d0/xrange
         do ipt=1,npt
            xs(ivar,ipt) = (xs(ivar,ipt) - xmin(ivar))*xrange
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------------
C Scale coordinates and function values
C-----------------------------------------------------------------------------
      subroutine scale_x(scalex, nvar, xes, xmin, xmax)
      implicit none
      include 'Param.h'
      integer          nvar, scalex
      double precision xes(nvar), xmin(nvar), xmax(nvar)

      integer          ivar
      double precision xrange
      
c If coordinate values dont need to be scaled
      if(scalex .eq. 0) return

c Find bounding box in design space and scale
      do ivar=1,nvar
         xrange = xmax(ivar) - xmin(ivar)
         if(xrange .lt. MACHEPS) then
            xrange = 1.0d0
         endif
         xrange = 1.0d0/xrange
         xes(ivar) = (xes(ivar) - xmin(ivar))*xrange
      enddo

      return
      end

C-----------------------------------------------------------------------------
C Scale coordinates and function values
C-----------------------------------------------------------------------------
      subroutine scale_fdata(scalef, npt, fs, fmin, fmax)
      implicit none
      include 'Param.h'
      integer          npt, scalef
      double precision fs(npt), fmin, fmax

      integer          ipt
      double precision frange
      
C     Find minimum and maximum function values
      fmin = +1.0d20
      fmax = -1.0d20
      do ipt=1,npt
         if(fs(ipt) .lt. fmin) fmin = fs(ipt)
         if(fs(ipt) .gt. fmax) fmax = fs(ipt)
      enddo

c If function values dont need to be scaled
      if(scalef .eq. 0)then
         print*,'Function will not be scaled'
         return
      endif

c Find range of function values and scale
      print*,'Scaling function values ...'
      frange = fmax - fmin
      if(frange .lt. MACHEPS)then
         write(*,*) 'Function is constant'
         frange = 1.0d0
      endif
      frange = 1.0d0/frange
      do ipt=1,npt
         fs(ipt) = (fs(ipt) - fmin)*frange
      enddo

      return
      end

C-----------------------------------------------------------------------------
C Unscale the data
C-----------------------------------------------------------------------------
      subroutine unscale_fun(scalef, fmin, fmax, fes)
      implicit none
      include 'Param.h'
      integer          scalef
      double precision fes, fmin, fmax

      double precision frange
      
c     If function values dont need to be scaled
      if(scalef .eq. 0) return

c     Find range of function values and scale
      frange = fmax - fmin
      if(frange .lt. MACHEPS)then
         frange = 1.0d0
      endif
      fes = fmin + frange*fes

      return
      end

C-----------------------------------------------------------------------------
C Unscale the data
C-----------------------------------------------------------------------------
      subroutine unscale_grad(scalex, scalef, nvar, xmin, xmax, 
     +                        fmin, fmax, fd)
      implicit none
      include 'Param.h'
      integer          scalex, scalef, nvar
      double precision fd(nvar), xmin(nvar), xmax(nvar), fmin, fmax

      integer          ivar
      double precision xrange, frange
      
      if(scalex.eq.0) goto 10

c     Scale coordinates
      do ivar=1,nvar
         xrange = xmax(ivar) - xmin(ivar)
         if(xrange .lt. MACHEPS) xrange = 1.0d0
         fd(ivar) = fd(ivar)/xrange
      enddo

10    if(scalef.eq.0) return

c     Scale function
      frange = fmax - fmin
      if(frange.lt.MACHEPS) frange = 1.0d0
      do ivar=1,nvar
         fd(ivar) = frange*fd(ivar)
      enddo

      return
      end

C-----------------------------------------------------------------------------
C Unscale the hessian
C-----------------------------------------------------------------------------
      subroutine unscale_hess(scalex, scalef, nvar, xmin, xmax, 
     +                        fmin, fmax, fdd)
      implicit none
      include 'Param.h'
      integer          scalex, scalef, nvar
      double precision fdd(nvar,nvar), xmin(nvar), xmax(nvar), 
     +                 fmin, fmax

      integer          ivar, jvar
      double precision xrange(nvar), fact, frange
      
      if(scalex.eq.0) goto 10

c     Scale coordinates
      do ivar=1,nvar
         xrange(ivar) = xmax(ivar) - xmin(ivar)
         if(xrange(ivar) .lt. MACHEPS) xrange(ivar) = 1.0d0
      enddo

      do ivar=1,nvar
         do jvar=1,nvar
            fact           = 1.0d0/(xrange(ivar)*xrange(jvar))
            fdd(ivar,jvar) = fact*fdd(ivar,jvar)
         enddo
      enddo

10    if(scalef.eq.0) return

c     Scale function value
      frange = fmax - fmin
      if(frange.lt.MACHEPS) frange = 1.0d0
      do ivar=1,nvar
         do jvar=1,nvar
            fdd(ivar,jvar) = frange*fdd(ivar,jvar)
         enddo
      enddo


      return
      end
