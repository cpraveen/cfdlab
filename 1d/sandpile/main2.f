      program main
      implicit none
      include 'vars.inc'

      real    ark(3)
      data    ark/0.0,0.75,0.33333333333333333333333/
      integer maxiter, nrk
      integer e, iter, rk
      real    xmin,xmax,dx,dt,umin,vmin,uxmax,uxc

      xmin    = 0.0
      xmax    = 1.0
      dx      = (xmax-xmin)/nel
      nrk     = 3
      maxiter = 1000000
      cfl     = 0.8

c     set initial condition
      do e=1,nel
         u(1,e) = 0.0
         u(2,e) = 0.0
         ua(e)  = 0.5*(u(1,e) + u(2,e))
         v(1,e) = 0.0
         v(2,e) = 0.0
         x(1,e) = xmin + (e-1)*dx
         x(2,e) = x(1,e) + dx
         B(e)   = x(1,e)
         B(e+1) = x(2,e)
      enddo
      call compute_ux()

      iter = 0
      do while(iter.lt.maxiter)
         uaold = ua
         uxold = ux
         vold  = v
         call compute_dt(dt)

         do rk=1,nrk
            call residual()

            ua = ark(rk)*uaold + (1.0-ark(rk))*(ua - dt*ures(1,:))
            ux = ark(rk)*uxold + (1.0-ark(rk))*(ux - dt*ures(2,:))
            v  = ark(rk)*vold  + (1.0-ark(rk))*(v  - dt*vres)

c           do e=1,nel
c              u(1,e) = ua(e) - 0.5*dx*ux(e)
c              u(2,e) = ua(e) + 0.5*dx*ux(e)
c           enddo
c           u(1,1)   = 0.0
c           u(2,nel) = 0.0

            u(1,1)   = 0.0
            u(2,1)   = dx*ux(1)
            do e=2,nel
               u(1,e) = u(2,e-1)
               u(2,e) = u(1,e) + dx*ux(e)
            enddo
            u(2,nel) = 0.0

            call limiter()
         enddo
         call check(umin, vmin, uxmax, uxc)

         iter = iter + 1
         write(*,10)iter,dt,umin,vmin,uxmax,uxc
10       format(i8,5e12.4)

         if(mod(iter,1000).eq.0) call output()
      enddo


      stop
      end
c------------------------------------------------------------------------------
      subroutine compute_dt(dt)
      implicit none
      include 'vars.inc'

      real    dt

      integer e
      real    dx

      dt = 1.0e20
      do e=1,nel
         dx = x(2,e) - x(1,e)
         dt = min(dt, dx/(v(1,e) + v(2,e) + 1.0e-20))
      enddo
      dt = cfl*dt
      dt = min(dt, 0.001)

      return
      end
c------------------------------------------------------------------------------
      subroutine compute_ux()
      implicit none
      include 'vars.inc'

      integer e

      do e=1,nel
         ux(e) = (u(2,e) - u(1,e))/(x(2,e) - x(1,e))
      enddo

      return
      end
c------------------------------------------------------------------------------
      subroutine num_flux_u(uxl, vl, uxr, vr, f)
      implicit none

      real uxl, uxr, vl, vr, f
      real f1, f2

      f1 = (abs(max(uxl,0.0)) - 1.0) * vl
      f2 = (abs(min(uxr,0.0)) - 1.0) * vr
      f  = max(f1, f2)

      return
      end
c------------------------------------------------------------------------------
      subroutine residual()
      implicit none
      include 'vars.inc'

      integer e
      real    f, vavg, Bavg, dx

c     u equation, only source term
      do e=1,nel
         vavg      = 0.5*(v(1,e) + v(2,e))
         ures(1,e) = -(1.0 - abs(ux(e)))*vavg
      enddo

      ures(2,:) = 0.0

c     flux at interior faces
      do e=1,nel-1
         call num_flux_u(ux(e), v(2,e), ux(e+1), v(1,e+1), f)
         ures(2,e)   = ures(2,e)   + f
         ures(2,e+1) = ures(2,e+1) - f
      enddo

c     v equation
      vres = 0.0

c     flux at left boundary
      e = 1
      f = ux(e)*v(1,e) + B(e)
      vres(1,e) = vres(1,e) + f

c     flux at interior faces
      do e=1,nel-1
         call num_flux(ux(e), v(2,e), ux(e+1), v(1,e+1), B(e+1), f)
         vres(2,e)   = vres(2,e)   - f
         vres(1,e+1) = vres(1,e+1) + f
      enddo

c     flux at right boundary
      e = nel
      f = ux(e)*v(2,e) + B(e+1)
      vres(2,e) = vres(2,e) - f

c     interior terms
      do e=1,nel
         vavg = 0.5*(v(1,e) + v(2,e))
         Bavg = 0.5*(B(e) + B(e+1))
         vres(1,e) = vres(1,e) - (vavg*ux(e)+Bavg)
         vres(2,e) = vres(2,e) + (vavg*ux(e)+Bavg)
      enddo

c     source term
      do e=1,nel
         dx        = x(2,e) - x(1,e)
         vres(1,e) = vres(1,e) + 0.5*dx*(1.0 - abs(ux(e)))*v(1,e)
         vres(2,e) = vres(2,e) + 0.5*dx*(1.0 - abs(ux(e)))*v(2,e)
         vres(:,e) = 2.0*vres(:,e)/dx
      enddo

      return
      end
c------------------------------------------------------------------------------
      subroutine num_flux(uxl, vl, uxr, vr, B, f)
      implicit none

      real uxl, uxr, vl, vr, B, f

      if(-uxl.ge.0.0.and.-uxr.gt.0.0)then
         f = uxl*vl + B
      else if(-uxl.lt.0.0.and.-uxr.le.0.0)then
         f = uxr*vr + B
      else if(-uxl.lt.0.0.and.-uxr.gt.0.0)then
         f = B
      else if(vl.gt.vr)then
         f = uxl*vl + B
      else if(vl.lt.vr)then
         f = uxr*vr + B
      else
         f = 0.5*(uxl*vl + uxr*vr) + B
      endif

      return
      end
c------------------------------------------------------------------------------
      real function minmod(a,b,c)
      implicit none
      real a, b, c

      if(a*b.gt.0.0.and.b*c.gt.0.0)then
         minmod = min(abs(a),abs(b))
         minmod = min(minmod,abs(c))
         if(a.le.0.0) minmod = -minmod
      else
         minmod = 0.0
      endif

      return
      end
c------------------------------------------------------------------------------
      subroutine limiter()
      implicit none
      include 'vars.inc'

      integer e
      real    ve, vem1, vep1, vx
      real    beta
      real    minmod

      beta = 1.0

      e = 1
      ve   = 0.5*(v(1,e  ) + v(2,e  ))
      vep1 = 0.5*(v(1,e+1) + v(2,e+1))
      vx   = v(2,e) - v(1,e)
      vx   = minmod(vx, beta*(vep1-ve), beta*(vep1-ve))
      v(1,e) = ve - 0.5*vx
      v(2,e) = ve + 0.5*vx

      do e=2,nel-1
         vem1 = 0.5*(v(1,e-1) + v(2,e-1))
         ve   = 0.5*(v(1,e  ) + v(2,e  ))
         vep1 = 0.5*(v(1,e+1) + v(2,e+1))
         vx   = v(2,e) - v(1,e)
         vx   = minmod(vx, beta*(ve-vem1), beta*(vep1-ve))
         v(1,e) = ve - 0.5*vx
         v(2,e) = ve + 0.5*vx
      enddo

      e = nel
      vem1 = 0.5*(v(1,e-1) + v(2,e-1))
      ve   = 0.5*(v(1,e  ) + v(2,e  ))
      vx   = v(2,e) - v(1,e)
      vx   = minmod(vx, beta*(ve-vem1), beta*(ve-vem1))
      v(1,e) = ve - 0.5*vx
      v(2,e) = ve + 0.5*vx

      return
      end
c------------------------------------------------------------------------------
      subroutine check(umin, vmin, uxmax, uxc)
      implicit none
      include 'vars.inc'

      real    umin, vmin, uxmax, uxc

      integer e
      real    ul, ur, dx

      umin = 1.0e20
      vmin = 1.0e20
      uxmax= 0.0
      do e=1,nel
         umin = min(umin, u(1,e))
         umin = min(umin, u(2,e))
         vmin = min(vmin, v(1,e))
         vmin = min(vmin, v(2,e))
         uxmax= max(uxmax, abs(ux(e)))
      enddo

      uxc  = 0.0
      do e=1,nel-1
         ul = 0.5*(u(1,e) + u(2,e))
         ur = 0.5*(u(1,e+1) + u(2,e+1))
         dx = x(2,e) - x(1,e)
         uxc = max(uxc, abs(ur - ul)/dx)
      enddo

      return
      end
c------------------------------------------------------------------------------
      subroutine output()
      implicit none
      include 'vars.inc'

      integer e, fid

      fid = 10
      open(fid, file='sol.dat')

      do e=1,nel
         write(fid,*) x(1,e), u(1,e), v(1,e)
         write(fid,*) x(2,e), u(2,e), v(2,e)
         write(fid,*)
      enddo

      close(fid)

      return
      end
