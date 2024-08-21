c     Solver 1d gas dynamic in Lagrangian form, non-conservative,
c     does not conserve total energy
c     Unknowns are
c        q(1) = specific volume = 1/density
c        q(2) = velocity 
c        q(3) = internal energy
      program main
      implicit none
      integer nvar, ncell
      parameter(nvar=3, ncell=1500)
      real gam, dx, dt, t, tf, xmin, xmax, xs, cfl
      real qold(nvar,ncell), q(nvar,ncell), res(nvar,ncell), x(ncell)
      real Dm(nvar), Dp(nvar), sl, sr, smax
      real vl, vr, ul, ur, pl, pr, p
      integer i, test, sod, shock, rk, nrk
      parameter(sod=1,shock=2,nrk=3)
      real ark(nrk), brk(nrk)
      ark = (/ 0.0, 3.0/4.0, 1.0/3.0 /)
      brk = (/ 1.0, 1.0/4.0, 2.0/3.0 /)

      gam  = 1.4
      xmin = 0.0
      xmax = 1.0
      xs   = 0.5
      cfl  = 0.1

      test = shock

      if(test.eq.sod)then
         vl = 1.0/1.000
         vr = 1.0/0.125
         ul = 0.0
         ur = 0.0
         pl = 1.0
         pr = 0.1
         tf = 0.2
      else if(test.eq.shock)then
c        Right moving shock: Abgrall/Karni, 2010
         vl = 1.0/0.4765625
         vr = 1.0/0.125
         ul = 2.304663838792128
         ur = 0.0
         pl = 1.0
         pr = 0.1
         tf = 0.5
      else
         stop "Unknown test"
      endif

      dx = (xmax - xmin)/ncell

c     Make grid, set ic
      do i=1,ncell
         x(i) = xmin + (i-1)*dx + 0.5*dx
         if(x(i) .lt. xs)then
            q(1,i) = vl
            q(2,i) = ul
            q(3,i) = pl*vl/(gam-1.0)
         else
            q(1,i) = vr
            q(2,i) = ur
            q(3,i) = pr*vr/(gam-1.0)
         endif
      enddo

c     Perform time stepping
      t = 0.0
      do while(t .lt. tf)
         qold = q
         do rk=1,nrk
            res = 0.0
            smax = 0.0
            do i=1,ncell-1 ! Loop only over interior faces
               call hll(gam, q(:,i), q(:,i+1), Dm, Dp, sl, sr)
               res(:,i  ) = res(:,i  ) + Dm
               res(:,i+1) = res(:,i+1) + Dp
               smax = max(smax,abs(sl))
               smax = max(smax,abs(sr))
            enddo
            if(rk.eq.1)then
               dt = cfl * dx / smax
               if(t + dt .gt. tf) dt = tf - t
            endif
            do i=1,ncell
               q(:,i) = ark(rk)*qold(:,i) 
     1                  + brk(rk)*(q(:,i) - (dt/dx)*res(:,i))
            enddo
         enddo
         t = t + dt
         print*,dt,t
      enddo

c     Save solution to file
      open(10,file="sol_noncon1.txt")
      write(10,*)"# Density, velocity, pressure"
      do i=1,ncell
         p = (gam-1.0)*q(3,i)/q(1,i)
         write(10,*)x(i),1.0/q(1,i),q(2,i),p
      enddo
      close(10)

      stop
      end

      subroutine hll(gam, ql, qr, Dm, Dp, sl, sr)
      implicit none
      integer nvar
      parameter(nvar=3)
      real gam, ql(nvar), qr(nvar), Dm(nvar), Dp(nvar), sl, sr
      real vl, vr, ul, ur, el, er, pl, pr, cl, cr
      real vs, us, es, ps, De

      vl = ql(1)
      ul = ql(2)
      el = ql(3)
      pl = (gam-1.0)*el/vl

      vr = qr(1)
      ur = qr(2)
      er = qr(3)
      pr = (gam-1.0)*er/vr

      cl = sqrt(gam * pl / vl)
      cr = sqrt(gam * pr / vr)
      sl = min(-cl, -cr)
      sr = max(cl, cr)

      vs = (sr*vr - sl*vl + (ur - ul))/(sr - sl)
      us = (sr*ur - sl*ul - (pr - pl))/(sr - sl)
      ps = pl + sl*(us - ul)

      Dm(1) = min(sl,0.0) * (vs - vl) + min(sr,0.0) * (vr - vs)
      Dm(2) = min(sl,0.0) * (us - ul) + min(sr,0.0) * (ur - us)

      Dp(1) = max(sl,0.0) * (vs - vl) + max(sr,0.0) * (vr - vs)
      Dp(2) = max(sl,0.0) * (us - ul) + max(sr,0.0) * (ur - us)

      De = 0.5*(pl + ps)*(us - ul) + 0.5*(ps + pr)*(ur - us)
      es = (sr*er - sl*el - De)/(sr - sl)
      Dm(3) = min(sl,0.0) * (es - el) + min(sr,0.0) * (er - es)
      Dp(3) = max(sl,0.0) * (es - el) + max(sr,0.0) * (er - es)

      return
      end subroutine hll
