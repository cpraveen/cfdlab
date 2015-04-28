c--------------------------------------------------------------------
c Solve linear convection equation with periodic bc
c Author: Praveen. C, http://praveen.tifrbng.res.in
c--------------------------------------------------------------------
      program main
      implicit none
      include 'common.inc'
      integer nc
      parameter(nc=100)
      real    f1(nc), f1old(nc), res1(nc)
      real    xmin, xmax, dx, x(nc), tf, cfl, dt
      integer i,irk, iter
      real    speed, maxspeed, t, f1min, f1max, f1tot

      xmin = 0.0  ! Left end of domain
      xmax = 1.0  ! Right end of domain
      tf   = 5.0  ! Final time

c     Set to ifirst/itvd/iweno3/iweno5
      recon_scheme = iweno5

c     cfl number
      if(recon_scheme.eq.ifirst) cfl = 0.9
      if(recon_scheme.eq.itvd) cfl = 0.45
      if(recon_scheme.eq.iweno3) cfl = 1.0/6.0
      if(recon_scheme.eq.iweno5) cfl = 1.0/12.0
      w1 = cfl

      dx   = (xmax - xmin)/nc

      maxspeed = 0.0
      do i=1,nc
         x(i) = 0.5*dx + (i-1)*dx
         maxspeed = max(maxspeed, abs(speed(x(i)-0.5*dx)))
      enddo
      dt = cfl * dx / maxspeed

      print*,"Maximum speed =", maxspeed
      print*,"Time Step     =", dt

c     set initial condition
      call set_initial_condition(nc, x, f1)
      call compute_range(nc, f1, dx, f1min, f1max, f1tot)
      print*,"Initial range =", f1min, f1max
      call output(nc, x, f1)
      call system("cp sol.dat init.dat")

      t = 0.0
      iter = 0
      do while(t.lt.tf)
         if(t+dt.gt.tf) dt = tf - t
         f1old = f1
         do irk=1,3
            call compute_residual (nc, dx, x, f1, res1)
            call update_solution (nc, f1, f1old, res1, irk, dt)
         enddo
         call compute_range(nc, f1, dx, f1min, f1max, f1tot)
         t = t + dt
         iter = iter + 1
         if(mod(iter,10).eq.0) call output(nc, x, f1)
         write(*,'(i8,4e14.6)') iter,t,f1min,f1max,f1tot
      enddo
      call output(nc, x, f1)

      stop
      end
c---------------------------------------------------------------------
c     TODO: We have to average initial condition
c---------------------------------------------------------------------
      subroutine set_initial_condition(nc, x, f1)
      implicit none
      integer nc
      real    x(nc), f1(nc)

      integer i
      real    pi
      parameter(pi=4.0*atan(1.0))

      do i=1,nc
         f1(i) = sin(2.0*pi*x(i))**2
      enddo

      return
      end
c---------------------------------------------------------------------
c Compute residual in the equation
c---------------------------------------------------------------------
      subroutine compute_residual (nc, dx, x, f1, res1)
      implicit none
      include 'common.inc'
      integer nc
      real    dx, x(nc), f1(nc), res1(nc)

      integer i
      real    f1l(nc+1), f1r(nc+1), flux1
      real    pstar, mini, maxi, theta1, theta2, theta

      res1(:) = 0.0

c     Compute interface values

c     First face
      call reconstruct(f1(nc-2), f1(nc-1), f1(nc), f1(1), f1(2), f1l(1))
      call reconstruct(f1(3),  f1(2),  f1(1),  f1(nc), f1(nc-1), f1r(1))

c     Second face
      call reconstruct(f1(nc-1), f1(nc), f1(1), f1(2), f1(3), f1l(2))
      call reconstruct(f1(4),  f1(3),  f1(2),  f1(1), f1(nc), f1r(2))

c     Third face
      call reconstruct(f1(nc), f1(1), f1(2), f1(3), f1(4), f1l(3))
      call reconstruct(f1(5),  f1(4),  f1(3),  f1(2), f1(1), f1r(3))

      do i=3,nc-3
         call reconstruct(f1(i-2), f1(i-1), f1(i), f1(i+1),
     1                    f1(i+2),f1l(i+1))
         call reconstruct(f1(i+3), f1(i+2), f1(i+1), f1(i),
     1                    f1(i-1),f1r(i+1))
      enddo

c     Third face from end
      call reconstruct(f1(nc-4),f1(nc-3),f1(nc-2),f1(nc-1),f1(nc),
     1                 f1l(nc-1))
      call reconstruct(f1(1), f1(nc), f1(nc-1), f1(nc-2),
     1                 f1(nc-3),f1r(nc-1))

c     penultimate face
      call reconstruct(f1(nc-3), f1(nc-2), f1(nc-1), f1(nc),
     1                 f1(1),f1l(nc))
      call reconstruct(f1(2), f1(1), f1(nc), f1(nc-1), f1(nc-2),f1r(nc))

c     Last face
      f1l(nc+1) = f1l(1)
      f1r(nc+1) = f1r(1)

c     Apply positivity limiter
      if(recon_scheme.eq.iweno3.or.recon_scheme.eq.iweno5)then
         do i=1,nc
            pstar = (f1(i) - w1*f1r(i) - w1*f1l(i+1))/(1.0 - 2*w1)
            mini = min(pstar, min(f1r(i), f1l(i+1)))
            maxi = max(pstar, max(f1r(i), f1l(i+1)))
            theta1 = abs(1.0-f1(i))/(abs(maxi-f1(i))+1.0e-14)
            theta2 = abs(0.0-f1(i))/(abs(mini-f1(i))+1.0e-14)
            theta  = min(1.0, min(theta1, theta2))
            f1r(i) = theta*(f1r(i) - f1(i)) + f1(i)
            f1l(i+1) = theta*(f1l(i+1) - f1(i)) + f1(i)
         enddo
      endif

c     Left state of first face and right state of last face
      f1l(1) = f1l(nc+1)
      f1r(nc+1) = f1r(1)

c     Compute fluxes

c     First face
      call num_flux(x(1)-0.5*dx, f1l(1), f1r(1), flux1)
      res1(nc) = res1(nc) + flux1
      res1(1)  = res1(1)  - flux1

      do i=1,nc-1
         call num_flux(x(i+1)-0.5*dx, f1l(i+1), f1r(i+1), flux1)
         res1(i)   = res1(i)   + flux1
         res1(i+1) = res1(i+1) - flux1
      enddo

      res1 = res1/dx

      return
      end
c---------------------------------------------------------------------
c Third order weno
c---------------------------------------------------------------------
      subroutine weno3(um1,u0,up1,u)
      implicit none
      real um1, u0, up1, u
      real eps, gamma1, gamma2
      parameter(eps = 1.0e-6, gamma1=1.0/3.0, gamma2=2.0/3.0)
      real beta1, beta2
      real u1, u2;
      real w1, w2;

      beta1 = (um1 - u0)**2
      beta2 = (up1 - u0)**2

      w1 = gamma1 / (eps+beta1)**2
      w2 = gamma2 / (eps+beta2)**2

      u1 = (3.0/2.0)*u0 - (1.0/2.0)*um1
      u2 = (u0 + up1)/2.0

      u = (w1 * u1 + w2 * u2)/(w1 + w2)

      return
      end
c---------------------------------------------------------------------
c Weno reconstruction
c---------------------------------------------------------------------
      subroutine weno5(um2,um1,u0,up1,up2,u)
      implicit none
      real um2, um1, u0, up1, up2, u
      real eps, gamma1, gamma2, gamma3
      parameter(eps = 1.0e-6, gamma1=1.0/10.0, gamma2=3.0/5.0,
     1          gamma3=3.0/10.0)
      real beta1, beta2, beta3
      real u1, u2, u3;
      real w1, w2, w3;

      beta1 = (13.0/12.0)*(um2 - 2.0*um1 + u0)**2 +
     1  (1.0/4.0)*(um2 - 4.0*um1 + 3.0*u0)**2
      beta2 = (13.0/12.0)*(um1 - 2.0*u0 + up1)**2 +
     1  (1.0/4.0)*(um1 - up1)**2
      beta3 = (13.0/12.0)*(u0 - 2.0*up1 + up2)**2 +
     1  (1.0/4.0)*(3.0*u0 - 4.0*up1 + up2)**2

      w1 = gamma1 / (eps+beta1)**2
      w2 = gamma2 / (eps+beta2)**2
      w3 = gamma3 / (eps+beta3)**2

      u1 = (1.0/3.0)*um2 - (7.0/6.0)*um1 + (11.0/6.0)*u0
      u2 = -(1.0/6.0)*um1 + (5.0/6.0)*u0 + (1.0/3.0)*up1
      u3 = (1.0/3.0)*u0 + (5.0/6.0)*up1 - (1.0/6.0)*up2

      u = (w1 * u1 + w2 * u2 + w3 * u3)/(w1 + w2 + w3)

      return
      end
c---------------------------------------------------------------------
      real function minmod(a, b, c)
      implicit none
      real a, b, c

      if(a*b.lt.0.0.or.b*c.lt.0.0)then
         minmod = 0.0
      else
         minmod = sign(1.0, a) * min(min(abs(a), abs(b)), abs(c))
      endif

      return
      end
c---------------------------------------------------------------------
      subroutine tvd(um2,um1,u0,up1,up2,u)
      implicit none
      real um2, um1, u0, up1, up2, u

      real dul, dur, duc, beta, minmod
      parameter(beta=2.0)

      dul = u0  - um1
      dur = up1 - u0
      duc = up1 - um1

      u = u0 + 0.5 * minmod(0.5*duc, beta*dul, beta*dur)

      return
      end
c---------------------------------------------------------------------
      subroutine reconstruct(um2,um1,u0,up1,up2,u)
      implicit none
      include 'common.inc'
      real um2, um1, u0, up1, up2, u

      if(recon_scheme.eq.ifirst)then
         u = u0
      else if(recon_scheme.eq.itvd)then
         call tvd(um2,um1,u0,up1,up2,u)
      else if(recon_scheme.eq.iweno3)then
         call weno3(um1,u0,up1,u)
      else if(recon_scheme.eq.iweno5)then
         call weno5(um2,um1,u0,up1,up2,u)
      else
         print*,"Unknown reconstruction scheme"
         stop
      endif

      return
      end
c---------------------------------------------------------------------
c Wave speed function
c---------------------------------------------------------------------
      real function speed (x)
      implicit none
      real x

      speed = 1.0

      return
      end
c---------------------------------------------------------------------
c Upwind Numerical flux
c---------------------------------------------------------------------
      subroutine num_flux(x, ul, ur, flux)
      implicit none
      real x, ul, ur, flux
      real speed, a

      a = speed(x)

      if(a.gt.0.0)then
         flux = a * ul
      else
         flux = a * ur
      endif

      return
      end
c---------------------------------------------------------------------
c 3-stage 3-order SSP-RK scheme
c---------------------------------------------------------------------
      subroutine update_solution (nc, f1, f1old, res1, irk, dt)
      implicit none
      integer nc, irk
      real    f1(nc), f1old(nc), res1(nc), dt

      if(irk.eq.1)then
         f1 = f1old - dt * res1
      else if(irk.eq.2)then
         f1 = (3.0/4.0)*f1old + (1.0/4.0)*(f1 - dt * res1)
      else if(irk.eq.3)then
         f1 = (1.0/3.0)*f1old + (2.0/3.0)*(f1 - dt * res1)
      endif

      return
      end
c---------------------------------------------------------------------
      subroutine compute_range(nc, f1, dx, f1min, f1max, f1tot)
      implicit none
      integer nc
      real    f1(nc), dx, f1min, f1max, f1tot

      integer i

      f1min = 1.0e20
      f1max =-1.0e20
      f1tot = 0.0

      do i=1,nc
         f1min = min(f1min, f1(i))
         f1max = max(f1max, f1(i))
         f1tot = f1tot + f1(i) * dx
      enddo

      return
      end
c---------------------------------------------------------------------
      subroutine output(nc, x, f1)
      implicit none
      integer nc
      real    x(nc), f1(nc)
      integer fid, i

      fid = 10
      open(fid,file="sol.dat")
      write(fid,'(2e24.12)')(x(i),f1(i),i=1,nc)
      close(fid)

      return
      end
