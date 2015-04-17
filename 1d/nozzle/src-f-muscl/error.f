      program main
      implicit none
      include 'param.inc'
      integer nf, nc, iter, i, j, nl, nr, pf
      integer nfmax, ncmax
      parameter(nfmax=1001)
      parameter(ncmax=nfmax-1)
      double precision x(nfmax), a(nfmax), q(3,ncmax), qold(3,ncmax), 
     +                 res(3,ncmax),
     +                 dt, residue, residue1, maxres, ar, cost, dum,
     +                 ptarg(ncmax)
      double precision qb(3,ncmax), resb(3,ncmax), resbold(3,ncmax), 
     +                 ab(nfmax), 
     +                 qb1(3,ncmax), costb
      integer          iflux1, iflux2, nc1, nc2, dnc
      double precision error1

      cfl = 0.8d0
      maxiter = 200000
      iflux1 = 1
      iflux2 = 1
      tolimit= 1
      noztyp = 1

      nc1 = 101
      nc2 = 101
      dnc = 100

      open(25,file='error.dat')
      do i=nc1,nc2,dnc

         open(20, file='param.dat')
         write(20,*) i
         write(20,*) 0.8
         write(20,*) maxiter
         write(20,*) iflux1, iflux2
         write(20,*) tolimit
         write(20,*) noztyp
         close(20)

         call system("../src-f-muscl/noz_flo > flo.log")
         call system("tail -n1 res.dat")

         open(20,file='cost.dat', status='old')
         read(20,*) cost
         close(20)


         call system("../src-f-muscl/noz_adj > flo.log")
         call system("tail -n1 resb.dat")

         call init(nf, nc)
         call read_flow(nc, q)
         call read_adj(nc, resb)
         call init_shape(nf, x, a)

c        call integrate1(nf, nc, q, resb, a, x, error1)
c        call integrate2(nf, nc, q, resb, a, x, error1)
         call integrate3(nf, nc, q, resb, a, x, error1)

         write(25,'(i6,2e24.15)') i, cost, cost+error1
         write(*,*) i, cost, cost+error1
         write(*,*)

      enddo
      close(25)

      stop
      end

      subroutine integrate1(nf, nc, q, resb, a, x, error1)
      implicit none
      include 'param.inc'
      integer nf, nc
      double precision q(3,*), resb(3,*), a(*), x(*), error1

      integer          i, k
      double precision x1, x2, h1, h2, u1, u2, p1, p2, e1, e2, e3
      double precision f1(3), f2(3)
      double precision y(3), w(3), g(3), z
      double precision xb
      double precision nozzle

c     nodes for 3-point gauss quadrature
      y(1) = -dsqrt(3.0d0/5.0d0)
      y(2) =  0.0d0
      y(3) = +dsqrt(3.0d0/5.0d0)

      w(1) = 5.0d0/9.0d0
      w(2) = 8.0d0/9.0d0
      w(3) = 5.0d0/9.0d0

      error1 = 0.0d0

      do i=2,nf-1

         x1 = 0.5d0*(x(i-1) + x(i))
         x2 = 0.5d0*(x(i+1) + x(i))

         h1 = nozzle(x1)
         h2 = nozzle(x2)

         u1 = q(2,i-1)/q(1,i-1)
         u2 = q(2,i  )/q(1,i  )

         p1 = gam1*( q(3,i-1) - 0.5d0*q(1,i-1)*u1**2 )
         p2 = gam1*( q(3,i  ) - 0.5d0*q(1,i  )*u2**2 )

         f1(1) = q(2,i-1)
         f2(1) = q(2,i  )

         f1(2) = p1 + q(1,i-1) * u1**2
         f2(2) = p2 + q(1,i  ) * u2**2

         f1(3) = (q(3,i-1) + p1)*u1
         f2(3) = (q(3,i  ) + p2)*u2

         e1 = 0.5d0*(resb(1,i-1) + resb(1,i)) * 
     1        ( h2*f2(1) - h1*f1(1) )

         e2 = 0.5d0*(resb(2,i-1) + resb(2,i)) *
     1        ( h2*f2(2) - h1*f1(2) )

         e3 = 0.5d0*(resb(3,i-1) + resb(3,i)) *
     1        ( h2*f2(3) - h1*f1(3) )

c        integrate v2(x)*h'(x)*p(x) using 3-point gauss
         do k=1,3

            z = 0.5d0*(x2-x1)*y(k) + 0.5d0*(x1+x2)

            call NOZZLE_X(L, ain, aout, z, xb, dc, dd, 1.0d0)
            g(k) = (resb(2,i-1)+(z-x1)*(resb(2,i)-resb(2,i-1))/(x2-x1))*
     1             (p1 + (z-x1)*(p2-p1)/(x2-x1)) * xb

         enddo

         e2 = e2 - 0.5d0*(x2-x1)*(w(1)*g(1) + w(2)*g(2) + w(3)*g(3))

         error1 = error1 + e1 + e2 + e3

      enddo

c        contribution from first cell
         p1 = gam1*(q(3,1) - 0.5d0*q(2,1)**2/q(1,1))
         z  = 0.5d0*(x(1) + x(2))
         h1 = nozzle(L, ain, aout, x(1), dc, dd)
         h2 = nozzle(L, ain, aout, z   , dc, dd)
         error1 = error1 - resb(2,1)*p1*(h2-h1)

c        contribution from last cell
         p1 = gam1*(q(3,nc) - 0.5d0*q(2,nc)**2/q(1,nc))
         z  = 0.5d0*(x(nf) + x(nf-1))
         h1 = nozzle(L, ain, aout, z   , dc, dd)
         h2 = nozzle(L, ain, aout, x(nf), dc, dd)
         error1 = error1 - resb(2,nc)*p1*(h2-h1)


      return
      end
c
c     error estimation using linear interpolation for Roe variables
c
      subroutine integrate2(nf, nc, q, resb, a, x, error1)
      implicit none
      include 'param.inc'
      integer nf, nc
      double precision q(3,*), resb(3,*), a(*), x(*), error1

      integer          i, j, k
      double precision z(3,nc+2), adj(3,nc+2), xg(nc+2), rho, u, p, H
      double precision x1, x2, xgauss, y(3), advar(3), w(3), errtmp
      double precision flux(3), fluxb(3), xgaussb, source
      double precision diverror, sourcerror

      do i=1,nc
         rho      = q(1,i)
         u        = q(2,i)/q(1,i)
         p        = gam1*( q(3,i) - 0.5d0*rho*u**2 )
         H        = gam*p/rho/gam1 + 0.5d0*u**2
         z(1,i+1) = dsqrt(rho)
         z(2,i+1) = dsqrt(rho) * u
         z(3,i+1) = dsqrt(rho) * H
         adj(:,i+1) = resb(:,i)
         xg(i+1) = 0.5d0*(x(i) + x(i+1))
      enddo

      xg(1) = x(1)
      xg(nc+2) = x(nc+1)

c     Extrapolate to two boundaries
      do i=1,3
         z(i,1   ) = 1.5d0*z(i,2   ) - 0.5d0*z(i,3 )
         z(i,nc+2) = 1.5d0*z(i,nc+1) - 0.5d0*z(i,nc)

         adj(i,1   ) = 1.5d0*adj(i,2   ) - 0.5d0*adj(i,3 )
         adj(i,nc+2) = 1.5d0*adj(i,nc+1) - 0.5d0*adj(i,nc)
      enddo

      open(22,file='z.dat')
      open(23,file='adj.dat')
      do i=1,nc+2
         write(22,'(4e16.8)')xg(i),z(1,i),z(2,i),z(3,i)
         write(23,'(4e16.8)')xg(i),adj(1,i),adj(2,i),adj(3,i)
      enddo
      close(22)
      close(23)

      y(1) = -dsqrt(3.0d0/5.0d0)
      y(2) = 0.0d0
      y(3) = +dsqrt(3.0d0/5.0d0)

      w(1) = 5.0d0/9.0d0
      w(2) = 8.0d0/9.0d0
      w(3) = 5.0d0/9.0d0

      open(30, file='errind.dat')
      error1 = 0.0d0
      do i=1,nc+1
         x1 = xg(i) 
         x2 = xg(i+1)
         errtmp = 0.0d0
         diverror=0.0d0
         sourcerror=0.0d0
         do j=1,3 ! loop over gauss points
            xgauss = 0.5d0*(x2 - x1)*y(j) + 0.5d0*(x1 + x2)
            do k=1,3
               advar(k) = adj(k,i) + 
     1            (xgauss-x1)*(adj(k,i+1)-adj(k,i))/(x2-x1)
            enddo
c           call InterpolateFlux(x1, x2, z(1,i), z(1,i+1), xgauss, flux)
            fluxb(:) = advar(:)
            call InterpolateFlux_x(x1, x2, z(1,i), z(1,i+1), 
     1                             xgauss, xgaussb, flux, fluxb)
            call InterpolateSource(x1, x2, z(1,i), z(1,i+1), 
     1                             xgauss, advar, source)
            diverror = diverror + w(j)*xgaussb
            sourcerror = sourcerror - w(j)*source
         enddo
         diverror = 0.5d0*(x2-x1)*diverror
         sourcerror = 0.5d0*(x2-x1)*sourcerror
         error1 = error1 + diverror + sourcerror
         write(30,'(4e16.8)') 0.5d0*(x1+x2), diverror, sourcerror,
     1               diverror+sourcerror
      enddo
      close(30)


      return
      end
c
c     error estimation using linear interpolation for Roe variables
c     first cell center at xmin, last at xmax
c     domain of integration is [xmin, xmax]
c
      subroutine integrate3(nf, nc, q, resb, a, x, error1)
      implicit none
      include 'param.inc'
      integer nf, nc
      double precision q(3,*), resb(3,*), a(*), x(*), error1

      integer          i, j, k
      double precision z(3,nc), adj(3,nc), xg(nc), rho, u, p, H
      double precision x1, x2, xgauss, y(3), advar(3), w(3), errtmp
      double precision flux(3), fluxb(3), xgaussb, source
      double precision diverror, sourcerror

      do i=1,nc
         rho      = q(1,i)
         u        = q(2,i)/q(1,i)
         p        = gam1*( q(3,i) - 0.5d0*rho*u**2 )
         H        = gam*p/rho/gam1 + 0.5d0*u**2
         z(1,i) = dsqrt(rho)
         z(2,i) = dsqrt(rho) * u
         z(3,i) = dsqrt(rho) * H
         adj(:,i) = resb(:,i)
         xg(i) = 0.5d0*(x(i) + x(i+1))
      enddo

      open(22,file='z.dat')
      open(23,file='adj.dat')
      do i=1,nc
         write(22,'(4e16.8)')xg(i),z(1,i),z(2,i),z(3,i)
         write(23,'(4e16.8)')xg(i),adj(1,i),adj(2,i),adj(3,i)
      enddo
      close(22)
      close(23)

      y(1) = -dsqrt(3.0d0/5.0d0)
      y(2) = 0.0d0
      y(3) = +dsqrt(3.0d0/5.0d0)

      w(1) = 5.0d0/9.0d0
      w(2) = 8.0d0/9.0d0
      w(3) = 5.0d0/9.0d0

      open(30, file='errind.dat')
      error1 = 0.0d0
      do i=1,nc-1
         x1 = xg(i) 
         x2 = xg(i+1)
         errtmp = 0.0d0
         diverror=0.0d0
         sourcerror=0.0d0
         do j=1,3 ! loop over gauss points
            xgauss = 0.5d0*(x2 - x1)*y(j) + 0.5d0*(x1 + x2)
            do k=1,3
               advar(k) = adj(k,i) + 
     1            (xgauss-x1)*(adj(k,i+1)-adj(k,i))/(x2-x1)
            enddo
c           call InterpolateFlux(x1, x2, z(1,i), z(1,i+1), xgauss, flux)
            fluxb(:) = advar(:)
            call InterpolateFlux_x(x1, x2, z(1,i), z(1,i+1), 
     1                             xgauss, xgaussb, flux, fluxb)
            call InterpolateSource(x1, x2, z(1,i), z(1,i+1), 
     1                             xgauss, advar, source)
            diverror = diverror + w(j)*xgaussb
            sourcerror = sourcerror - w(j)*source
         enddo
         diverror = 0.5d0*(x2-x1)*diverror
         sourcerror = 0.5d0*(x2-x1)*sourcerror
         error1 = error1 + diverror + sourcerror
         write(30,'(4e16.8)') 0.5d0*(x1+x2), diverror, sourcerror,
     1               diverror+sourcerror
      enddo
      close(30)


      return
      end
