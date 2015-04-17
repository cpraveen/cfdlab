      program main
      implicit none
      integer          nf, i
      parameter(nf=101)
      real x(nf), a(nf), ab(nf), adum(nf), grad

      open(10, file='shape.dat')
      do i=1,nf
         read(10,*) x(i), a(i)
         adum(i) = ab(i)
      enddo
      close(10)

      open(10, file='shapeb.dat')
      do i=1,nf
         read(10,*) x(i), ab(i)
         adum(i) = ab(i)
      enddo
      close(10)

      grad = 0.0
      do i=1,nf
         grad = grad + ab(i)**2
      enddo
      grad = sqrt(grad)
      print*,'Shape gradient =',grad
c     do i=1,nf
c        ab(i) = ab(i)/grad
c     enddo
      call smooth_shape(nf, x, ab)

      open(10, file='out')
      do i=1,nf
         write(10,*) x(i), a(i) - 0.05*ab(i), ab(i)
      enddo
      close(10)

      stop
      end

C Subroutine to smooth the surface deformation by solving an elliptic
C equation following Mohammadi and Pironneau.
C Single precision version
      subroutine smooth_shape(n, x, dx)
      implicit none
      integer n
      real    x(n), dx(n)

      integer i, iter
      real    residue, residue1, RTOL, s(n), 
     &        dx1(n), dxo(n), eps(n)
      parameter(RTOL=1.0e-6)

      print*,'Smoothening shape ...'

c First point is assumed to be trailing edge which is fixed
c These things must be specified in an input file in future
      dx(1) = 0.0
      dx(n) = 0.0

      residue = 1.0
      iter    = 0

      call saveshape(n, dx, dxo)
      call epsi(n, dxo, eps)
      open(unit=10, file='ress.dat')
      do while(residue .gt. RTOL .and. iter .lt. 1000)
         iter = iter + 1
         call saveshape(n, dx, dx1)
         call smoothiter2(n, dx, dxo, dx1, x, eps)
         call smoothres(n, dx, dx1, iter, residue, residue1)
         write(10,*) iter,residue
      enddo
      close(10)

      return
      end

C Viscosity coefficient for smoothing. Viscosity is set to zero if the
C local shape is smooth
      subroutine epsi(n, dx, eps)
      implicit none
      integer n
      real    dx(n), eps(n)

      integer i
      real    dxa, f, dl, dr, dc, ELIMIT, am

c Find average displacement of shape, using L1 norm
      dxa = 0.0
      do i=1,n
c        dxa = dxa + abs(dx(i))
c        dxa = dxa + abs(dx(i))**2
         dxa = amax1(dxa, abs(dx(i)))
      enddo
c     dxa = dxa/n
c     dxa = sqrt(dxa/n)

      dxa = 0.0
      do i=2,n-1
         dl  = dx(i) - dx(i-1)
         dr  = dx(i+1) - dx(i)
         dxa = amax1(dxa, abs(dl)+abs(dr))
      enddo

      print*,'Average shape change:'
      print*,'\t dx =',dxa

      if(dxa .lt. 1.0e-10) stop

      do i=1,n
         eps(i) = 0.0
      enddo

      do i=2,n-1
         dl = dx(i)   - dx(i-1)
         dr = dx(i+1) - dx(i)
         dc = dx(i+1) - dx(i-1)

         f = (abs(dr) + abs(dl))
         eps(i) = 1.0*f/dxa

c        f = abs( abs(dl) - abs(dr) )
c        am= amax1( abs(dl)+abs(dr), 1.0e-10 )
c        eps(i) = f/am

         print*,i,f,eps(i)
      enddo

      return
      end

C Save the current shape change into another array
      subroutine saveshape(n, dx, dxo)
      implicit none
      integer n
      real    dx(n), dxo(n)

      integer i

      do i=1,n
         dxo(i) = dx(i)
      enddo

      return
      end

C Perform one iteration of smoothing
C Finite difference scheme
      subroutine smoothiter2(n, dx, dxo, dx1, s, eps)
      implicit none
      integer n
      real    dx(n), dxo(n), dx1(n), s(n), eps(n)

      integer i
      real    h, epsl, epsr, mo

      do i=2,n-1
         epsl  = eps(i)
         epsr  = eps(i)
         h     = 0.5*( s(i+1) - s(i-1) )
         mo    = h + epsl/(s(i)-s(i-1)) + epsr/(s(i+1)-s(i))
         dx(i) = h*dxo(i) + epsl*dx1(i-1)/(s(i)-s(i-1)) +
     &                      epsr*dx1(i+1)/(s(i+1)-s(i))
         dx(i) = dx(i)/mo
      enddo

      return
      end

C Perform one iteration of smoothing
C Finite volume scheme
      subroutine smoothiter3(n, dx, dxo, dx1, s, eps)
      implicit none
      integer n
      real    dx(n), dxo(n), dx1(n), s(n), eps(n)

      integer i
      real    h, epsl, epsr, mo

      do i=2,n-1
c        epsl  = 0.5*( eps(i-1) + eps(i) )
c        epsr  = 0.5*( eps(i+1) + eps(i) )
         epsl  = amax1( eps(i-1) , eps(i) )
         epsr  = amax1( eps(i+1) , eps(i) )
         h     = 0.5*( s(i+1) - s(i-1) )
         mo    = h + epsl/(s(i)-s(i-1)) + epsr/(s(i+1)-s(i))
         dx(i) = h*dxo(i) + epsl*dx1(i-1)/(s(i)-s(i-1)) +
     &                      epsr*dx1(i+1)/(s(i+1)-s(i))
         dx(i) = dx(i)/mo
      enddo

      return
      end

C Find residue for change after iteration 
      subroutine smoothres(n, dx, dx1, iter, residue, residue1)
      implicit none
      integer n, iter
      real    dx(n), dx1(n), residue, residue1

      integer i

      residue = 0.0
      do i=1,n
         residue = residue + (dx(i) - dx1(i))**2
      enddo
      residue = sqrt(residue)/n

      if(iter .eq. 1)then
         residue1 = residue
         print*,'Residue in first iteration =',residue1
      endif

      residue = residue/residue1

      return
      end

      real function ELIMIT(dl, dr, dc)
      implicit none
      real dl, dr, dc, top, bot, EPSILON
      EPSILON = 1.0e-15

      top = dr - dl
      bot = amax1( abs(dl)+abs(dr), 1.0e-2 )
      if( dl*dr .le. 0.0)then
         elimit = 0.0
      else
         elimit = 1.0
      endif

      return
      end
