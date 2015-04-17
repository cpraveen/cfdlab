C Initialize some constants and grid
      subroutine init(nf, nc)
      implicit none
      include 'param.inc'
      integer nf, nc

      ark(1) = 0.0d0
      ark(2) = 1.0/3.0d0
      ark(3) = 2.0d0/3.0d0

      gam  = 1.4d0
      gam1 = gam - 1.0d0
      gasconst = 1.0d0

      minf = 1.5d0
      rinf = 1.0d0
      pinf = 1.0d0
      Tinf = pinf/(rinf*gasconst)
      uinf = minf*dsqrt(gam*gasconst*Tinf)

      prat = 2.5d0
      pout = prat*pinf

      qinf(1) = rinf
      qinf(2) = rinf*uinf
      qinf(3) = pinf/gam1 + 0.5d0*rinf*uinf**2

      print*,'Inflow mach number =',minf
      print*,'Pressure ratio     =',prat

      print*,'Reading paramters from param.dat'
      open(50, file="param.dat")
      read(50,*) nc
      read(50,*) cfl
      read(50,*) maxiter
      read(50,*) flux1, flux2
      read(50,*) tolimit
      read(50,*) noztyp
      close(50)

      nf = nc + 1

      return
      end

C Copy array q into qold
      subroutine save_old(nc, q, qold)
      implicit none
      integer          nc
      double precision q(3,*), qold(3,*)

      integer          i, j

      do i=1,nc
         do j=1,3
            qold(j,i) = q(j,i)
         enddo
      enddo

      return
      end

C Convert conservative to primitive variables
      subroutine con_to_prim(r, u, p, q)
      implicit none
      include 'param.inc'
      double precision r, u, p, q(3)

      r = q(1)
      u = q(2)/q(1)
      p = gam1*( q(3) - 0.5d0*q(2)**2/q(1) )

      return
      end

C Global time step
      subroutine timestep(nc, q, dt)
      implicit none
      include 'param.inc'
      integer          nc
      double precision q(3,*), dt

      integer          i
      double precision dtlocal, a, r, u, p

      dt = 1.0d20
      do i=1,nc
         call con_to_prim(r, u, p, q(1,i))
         a = dsqrt(gam*p/r)
         dtlocal = dx/( dabs(u) + a )
         dt = dmin1(dt, dtlocal)
      enddo
      dt = cfl*dt

      return
      end
     
C L2 norm of solution residual
      subroutine sol_residue(nc, res, residue)
      implicit none
      integer          nc
      double precision res(3,*), residue

      integer          i

      residue = 0.0d0
      do i=1,nc
         residue = residue + res(1,i)**2 + res(2,i)**2 + res(3,i)**2
      enddo
      residue = dsqrt(residue)/nc

      return
      end

C Save flow solution and nozzle shape
      subroutine result(nf, nc, x, q)
      implicit none
      include 'param.inc'
      integer          nf, nc
      double precision x(*), q(3,*)

      integer          i
      double precision r, u, p, son, mach

      open(unit=11, file='flow.dat')
      do i=1,nc
         call con_to_prim(r, u, p, q(1,i))
         son  = dsqrt(gam*p/r)
         mach = u/son
         write(11,'(5e24.15)') x(i)+0.5d0*dx, r, u, p, mach
      enddo
      close(11)

      return
      end

C Read flow solution
      subroutine read_flow(nc, q)
      implicit none
      include 'param.inc'
      integer          nc
      double precision q(3,*)

      integer          i
      double precision r, u, p, dum, mach

      print*,'Reading flow solution from flow.dat'
      open(unit=11, file='flow.dat', status='old')
      do i=1,nc
         read(11,*) dum, r, u, p, mach
         q(1,i) = r
         q(2,i) = r*u
         q(3,i) = p/gam1 + 0.5d0*r*u**2
      enddo
      close(11)

      return
      end

C Read adjoint solution
      subroutine read_adj(nc, q)
      implicit none
      include 'param.inc'
      integer          nc
      double precision q(3,*)

      integer          i
      double precision r, u, p, dum, mach

      print*,'Reading adjoint solution from flowb.dat'
      open(unit=11, file='flowb.dat', status='old')
      do i=1,nc
         read(11,*) dum, q(1,i), q(2,i), q(3,i)
      enddo
      close(11)

      return
      end

C Read nozzle shape
      subroutine read_shape(nf, x, a)
      implicit none
      include 'param.inc'
      integer          nf
      double precision x(*), a(*)

      integer          i

      print*,'Reading nozzle shape from shape.dat'
      open(unit=11, file='shape.dat')
      do i=1,nf
         read(11,*) x(i), a(i)
      enddo
      close(11)

      dx   = x(2) - x(1)
      ain  = a(1)
      aout = a(nf)

      print*,'Inlet area         =',ain
      print*,'Outlet area        =',aout
      print*,'Area ratio         =',aout/ain

      return
      end

C Save adjoint flow solution and shape gradient
      subroutine result_adj(nf, nc, x, a, q)
      implicit none
      include 'param.inc'
      integer          nf, nc
      double precision x(*), a(*), q(3,*)

      integer          i
      double precision r, u, p, son, mach

      open(unit=11, file='flowb.dat')
      do i=1,nc
         write(11,'(4e24.15)') x(i)+0.5d0*dx, q(1,i), q(2,i), q(3,i)
      enddo
      close(11)

      open(unit=11, file='shapeb.dat')
      do i=1,nf
         write(11,'(2e24.15)') x(i), a(i)
      enddo
      close(11)

      return
      end

C Set nozzle shape
      subroutine init_shape(nf, x, a)
      implicit none
      include 'param.inc'
      integer          nf
      double precision x(*), a(*)

      integer          i
      double precision nozzle

      if(noztyp.eq.1.or.noztyp.eq.3)then
         L   = 10.0d0
         ain = 1.0512d0
         aout= 1.75d0
         dc  = 0.8d0
         dd  = 4.0d0
         db  = (aout - ain)/(dtanh(10.d0*dc - dd) - dtanh(-dd))
         da  = ain - db*dtanh(-dd)
         xmin= 0.0d0
      elseif(noztyp.eq.2)then
         L   = 2.0d0
         ain = 2.0d0
         aout= 2.0d0
         xmin=-1.0d0
      elseif(noztyp.eq.4)then
         L   = 4.0d0
         xmin= -2.0d0
         ain = 2.0d0
         ath = 1.0d0
      else
         print*,'nozzle type unknown'
         stop
      endif

c     First cell face is at x=xmin - 0.5*dx
c     last face at x = xmax + 0.5*dx
c     first cell center is at x = xmin and last cell center at x = xmax
      dx  = L/(nf-2)
      open(unit=10, file='shape.dat')
      do i=1,nf
         x(i) = xmin + dx*(i-1) - 0.5*dx
         a(i) = nozzle(x(i))
         write(10,'(2e20.10)') x(i), a(i)
      enddo
      close(10)

      print*,'Number of points   =',nf
      print*,'Length of nozzle   =',L
      print*,'dx                 =',dx
      print*,'Inlet area         =',ain
      print*,'Outlet area        =',aout
      print*,'Area ratio         =',aout/ain
      print*,'a                  =',da
      print*,'b                  =',db
      print*,'c                  =',dc
      print*,'d                  =',dd

      print*,'Shape written into shape.dat'

      return
      end
