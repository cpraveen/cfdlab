

C Flux for inflow face
C Use inflow conditions to compute flux
      subroutine flux_in(a, res)
      implicit none
      include 'param.inc'
      double precision a, res(3)

      integer          i
      double precision einf, fl(3)

      fl(1) = rinf * uinf
      fl(2) = pinf + rinf * uinf * uinf
      einf  = pinf/gam1 + 0.5d0*rinf*uinf**2
      fl(3) = (einf + pinf) * uinf

      do i=1,3
         res(i) = res(i) - a*fl(i)
      enddo

      return
      end

C Flux for inflow face
C Use rho and p from inflow, velocity from inside
      subroutine flux_in2(a, q, res)
      implicit none
      include 'param.inc'
      double precision a, q(3), res(3)

      integer          i
      double precision einf, fl(3), u

      u     = q(2)/q(1)

      fl(1) = rinf * u
      fl(2) = pinf + rinf * u * u
c     einf  = pinf/gam1 + 0.5d0*rinf*u**2
      einf  = pinf/gam1 + 0.5d0*rinf*uinf**2
      fl(3) = (einf + pinf) * u

      do i=1,3
         res(i) = res(i) - a*fl(i)
      enddo

      return
      end

c     Inflow flux based on steger warming splitting
      subroutine flux_in_steger(a, q, res)
      implicit none
      include 'param.inc'
      double precision a, q(3), res(3)

      integer          i
      double precision rl, ul, pl, cl, hl, l1, l2, l3, f(3), g(3)
      double precision rr, ur, pr, cr, hr, r1, r2, r3

c     positive flux from left
      rl = rinf
      pl = pinf
      cl = dsqrt(gam*pl/rl)
      ul = minf*cl
      hl = cl**2/gam1 + 0.5d0*ul**2

      l1 = max(0.0d0, ul-cl)
      l2 = max(0.0d0, ul)
      l3 = max(0.0d0, ul+cl)

      f(1) = (0.5d0*rl/gam)*(l1 + 2.0d0*gam1*l2 + l3 )
      f(2) = (0.5d0*rl/gam)*((ul-cl)*l1 + 2.0d0*gam1*ul*l2 + (ul+cl)*l3)
      f(3) = (0.5d0*rl/gam)*((hl-ul*cl)*l1 + gam1*ul**2*l2 +
     1                       (hl+ul*cl)*l3)

c     negative flux from right side
      rr = q(1)
      ur = q(2)/q(1)
      pr = gam1*(q(3) - 0.5d0*q(2)**2/q(1))
      cr = dsqrt(gam*pr/rr)
      hr = cr**2/gam1 + 0.5d0*ur**2

      r1 = min(0.0d0, ur-cr)
      r2 = min(0.0d0, ur)
      r3 = min(0.0d0, ur+cr)

      g(1) = (0.5d0*rr/gam)*(r1 + 2.0d0*gam1*r2 + r3 )
      g(2) = (0.5d0*rr/gam)*((ur-cr)*r1 + 2.0d0*gam1*ur*r2 + (ur+cr)*r3)
      g(3) = (0.5d0*rr/gam)*((hr-ur*cr)*r1 + gam1*ur**2*r2 +
     1                       (hr+ur*cr)*r3)

      do i=1,3
         res(i) = res(i) - a*(f(i) + g(i))
      enddo

      return
      end

C Flux for interior face
C AUSMDV flux, taken from nozzle code of Manoj Nair, CTFD, NAL
      subroutine flux_ausm(a, ql, qr, resl, resr)
      implicit none
      include 'param.inc'
      double precision a, ql(3), qr(3), resl(3), resr(3)

      integer          i
      double precision rho_left, rho_right, u_left, u_right,
     +                 p_left, p_right, c_left, c_right, c_mean,
     +                 H_left, H_right, amach_l, amach_r, 
     +                 alpha_l, alpha_r, temp1, temp2, u_lp, u_rm,
     +                 p_lp, p_rm, ak_half, p_half, s_half, flux_v1,
     +                 flux_d1, flux_p, flux_n, fl(3)

      rho_left = ql(1)
      u_left   = ql(2)/ql(1)
      p_left   = gam1*( ql(3) - 0.5d0*ql(2)**2/ql(1) )

      rho_right = qr(1)
      u_right   = qr(2)/qr(1)
      p_right   = gam1*( qr(3) - 0.5d0*qr(2)**2/qr(1) )

      C_LEFT   = dsqrt(gam*p_left/rho_left)
      H_LEFT   = c_left**2/gam1 + 0.5d0*u_left**2
    
C.......Right State 
      C_RIGHT  = dsqrt(gam*p_right/rho_right)
      H_RIGHT  = c_right**2/gam1 + 0.5d0*u_right**2

C.......Mean Speed of Sound
      C_MEAN = dmax1(C_LEFT,C_RIGHT)

C.......Left and Right Mach no.
      AMACH_L = DABS(U_LEFT)/C_MEAN
      AMACH_R = DABS(U_RIGHT)/C_MEAN

C.......ALPHA
      TEMP1   = P_LEFT/RHO_LEFT
      TEMP2   = P_RIGHT/RHO_RIGHT

      ALPHA_L = 2*TEMP1/(TEMP1 + TEMP2)
      ALPHA_R = 2*TEMP2/(TEMP1 + TEMP2)

C.......U_LEFT+
      TEMP1 = U_LEFT + C_MEAN
      TEMP2 = U_LEFT + DABS(U_LEFT)
      IF(AMACH_L .LE. 1) THEN
        U_LP  = ALPHA_L*(0.25d0*TEMP1*TEMP1/C_MEAN - 0.5d0*TEMP2)
     &               + 0.5d0*TEMP2
       ELSE
        U_LP  = 0.5d0*TEMP2
      ENDIF
    
C.......U_RIGHT-
      TEMP1 = U_RIGHT - C_MEAN
      TEMP2 = U_RIGHT - DABS(U_RIGHT)
      IF(AMACH_R .LE. 1) THEN
        U_RM  = ALPHA_R*(-0.25d0*TEMP1*TEMP1/C_MEAN - 0.5d0*TEMP2)
     &               + 0.5d0*TEMP2
       ELSE
        U_RM  = 0.5d0*TEMP2
      ENDIF

C.......P_LEFT+
      TEMP1 = U_LEFT/C_MEAN + 1
      TEMP2 = 2 - U_LEFT/C_MEAN
      IF(AMACH_L .LE. 1) THEN
        P_LP = 0.25d0*P_LEFT*TEMP1*TEMP1*TEMP2
       ELSE
        P_LP = 0.5d0*P_LEFT*(U_LEFT + DABS(U_LEFT))/(U_LEFT)
      ENDIF

C.......P_RIGHT-
      TEMP1 = U_RIGHT/C_MEAN - 1
      TEMP2 = 2 + U_RIGHT/C_MEAN
      IF(AMACH_R .LE. 1) THEN
        P_RM = 0.25d0*P_RIGHT*TEMP1*TEMP1*TEMP2
       ELSE
        P_RM = 0.5d0*P_RIGHT*(U_RIGHT - DABS(U_RIGHT))/(U_RIGHT)
      ENDIF

      AK_HALF = 10.d0
      P_HALF  = P_LP + P_RM
      S_HALF  = DABS(P_RIGHT-P_LEFT)/dmin1(P_LEFT,P_RIGHT)
      S_HALF  = 0.5d0*dmin1(1.d0,AK_HALF*S_HALF)

C.......Calculating the inviscid interface fluxes

C.......Mass flux
      FL(1) = U_LP*RHO_LEFT + U_RM*RHO_RIGHT
    
C.......Enthalpy flux
      FL(3) = 0.5d0*(FL(1)*(H_LEFT + H_RIGHT) 
     &               - DABS(FL(1))*(H_RIGHT - H_LEFT))
 
C.......Momentum flux

      FLUX_V1 = RHO_LEFT*U_LEFT*U_LP + RHO_RIGHT*U_RIGHT*U_RM

      FLUX_D1 = 0.5d0*(FL(1)*(U_LEFT + U_RIGHT) - 
     &                DABS(FL(1))*(U_RIGHT - U_LEFT))

      FLUX_P = P_HALF

      FLUX_N  = (0.5d0+S_HALF)*FLUX_V1 + (0.5d0-S_HALF)*FLUX_D1

      FL(2) = FLUX_N  + FLUX_P

      do i=1,3
         resl(i) = resl(i) + a*fl(i)
         resr(i) = resr(i) - a*fl(i)
      enddo

      return
      end

C Flux for interior face
C Lax-Fredrichs flux
      subroutine flux_lf(a, ql, qr, resl, resr)
      implicit none
      include 'param.inc'
      double precision a, ql(3), qr(3), resl(3), resr(3)

      integer          i
      double precision rho_left, rho_right, u_left, u_right,
     +                 p_left, p_right, c_left, c_right,
     +                 H_left, H_right, rl, rr, rs, ua, ha, aa, lam,
     +                 fl(3), fr(3), flux(3)

      rho_left = ql(1)
      u_left   = ql(2)/ql(1)
      p_left   = gam1*( ql(3) - 0.5d0*ql(2)**2/ql(1) )
      C_LEFT   = dsqrt(gam*p_left/rho_left)
      H_LEFT   = c_left**2/gam1 + 0.5d0*u_left**2
    

C.......Right State 
      rho_right = qr(1)
      u_right   = qr(2)/qr(1)
      p_right   = gam1*( qr(3) - 0.5d0*qr(2)**2/qr(1) )
      C_RIGHT  = dsqrt(gam*p_right/rho_right)
      H_RIGHT  = c_right**2/gam1 + 0.5d0*u_right**2

      fl(1) = ql(2)
      fl(2) = p_left + rho_left*u_left**2
      fl(3) = ql(2)*H_left

      fr(1) = qr(2)
      fr(2) = p_right + rho_right*u_right**2
      fr(3) = qr(2)*H_right

C.......Roe averages
      rl = dsqrt(rho_left)
      rr = dsqrt(rho_right)
      rs = 1.0d0/(rl + rr)
      ua = (rl*u_left + rr*u_right)*rs
      ha = (rl*h_left + rr*h_right)*rs
      aa = dsqrt( (gam-1.0d0)*( ha - 0.5d0*ua**2 ) )

c     ua = 0.5d0*( u_left + u_right)
c     aa = 0.5d0*( c_left + c_right)

      lam= dabs(ua) + aa

      do i=1,3
         flux(i) = 0.5d0*( fl(i) + fr(i) ) - 0.5d0*lam*( qr(i) - ql(i) )
         resl(i) = resl(i) + a*flux(i)
         resr(i) = resr(i) - a*flux(i)
      enddo

      return
      end

C Flux for interior face
C KFVS flux, taken from nozzle code of Manoj Nair
      subroutine flux_kfvs(a, ql, qr, resl, resr)
      implicit none
      include 'param.inc'
      double precision a, ql(3), qr(3), resl(3), resr(3)

      integer          i
      double precision rl, rr, ul, ur, pl, pr, el, er, betal, betar,
     +                 sl, sr, Al, Ar, Bl, Br, Fp(3), Fm(3), ERRF

      rl = ql(1)
      ul = ql(2)/ql(1)
      pl = gam1*( ql(3) - 0.5d0*ql(2)**2/ql(1) )
      el = pl/GAM1 + 0.5d0*rl*ul**2

      rr = qr(1)
      ur = qr(2)/qr(1)
      pr = gam1*( qr(3) - 0.5d0*qr(2)**2/qr(1) )
      er = pr/GAM1 + 0.5d0*rr*ur**2

      betal = 0.5d0*rl/pl
      sl    = ul*dsqrt(betal)
      Al    = 0.5d0*(1.0d0 + ERRF(sl))
      Bl    = 0.5d0*dexp(-sl**2)/dsqrt(M_PI*betal)

      Fp(1) = rl*(ul*Al + Bl)
      Fp(2) = (pl + rl*ul**2)*Al + rl*ul*Bl
      Fp(3) = (el + pl)*ul*Al + (el + 0.5d0*pl)*Bl

c     Negative flux
      betar = 0.5d0*rr/pr
      sr    = ur*dsqrt(betar)
      Ar    = 0.5d0*(1.0d0 - ERRF(sr))
      Br    = 0.5d0*dexp(-sr**2)/dsqrt(M_PI*betar)

      Fm(1) = rr*(ur*Ar - Br)
      Fm(2) = (pr + rr*ur**2)*Ar - rr*ur*Br
      Fm(3) = (er + pr)*ur*Ar - (er + 0.5d0*pr)*Br

      do i=1,3
         resl(i) = resl(i) + a*(fp(i) + fm(i))
         resr(i) = resr(i) - a*(fp(i) + fm(i))
      enddo

      return
      end

C Flux for outflow face
C Linear extrapolation from last two cells to the outflow face
C Pressure taken from pout
      subroutine flux_out(a, ql, qll, res)
      implicit none
      include 'param.inc'
      double precision a, ql(3), qll(3), res(3)

      integer          i
      double precision rl, ul, rll, ull
      double precision r, u, p, e, fl(3)

      rl  = ql(1)
      ul  = ql(2)/ql(1)

      rll = qll(1)
      ull = qll(2)/qll(1)

      r = 1.5d0*rl - 0.5d0*rll
      u = 1.5d0*ul - 0.5d0*ull
      e = pout/gam1 + 0.5d0*r*u**2

      fl(1) = r*u
      fl(2) = pout + r*u**2
      fl(3) = (e + pout)*u

      do i=1,3
         res(i) = res(i) + a*fl(i)
      enddo

      return
      end

C Flux for outflow face
      subroutine flux_out2(a, ql, res)
      implicit none
      include 'param.inc'
      double precision a, ql(3), res(3)

      integer          i
      double precision r, u, p, e, fl(3)
      double precision rl, ul, pl
      double precision c
      double precision Sl, Hl, cl

      rl = ql(1)
      ul = ql(2)/ql(1)
      pl = gam1*(ql(3) - 0.5d0*rl*ul*ul)
      Sl = pl/rl**gam
      cl = dsqrt(gam*pl/rl)
      Hl = ul + 2.0d0*cl/gam1

      r  = (pout/Sl)**(1.0d0/gam)
      c  = dsqrt(gam*pout/r)
      u  = Hl - 2.0d0*c/gam1

      e = pout/gam1 + 0.5d0*r*u**2
      fl(1) = r*u
      fl(2) = pout + r*u**2
      fl(3) = (e + pout)*u
      
c     e = pl/gam1 + 0.5d0*rl*ul**2
c     fl(1) = rl*ul
c     fl(2) = pl + rl*ul**2
c     fl(3) = (e + pl)*ul

      do i=1,3
         res(i) = res(i) + a*fl(i)
      enddo

      return
      end

C Source term contribution to residual
      subroutine source_term(al, ar, q, res)
      implicit none
      include 'param.inc'
      double precision al, ar, q(3), res(3)

      double precision p

      p      = gam1*( q(3) - 0.5d0*q(2)**2/q(1) )
      res(2) = res(2) - p*( ar - al )

      return
      end

C Cost function, L2 norm of pressure difference
      subroutine costfunc(nc, q, cost)
      implicit none
      include 'param.inc'
      integer          nc
      double precision q(3,*), cost

      double precision p, p1, p2
      integer          i

      cost = 0.0d0
      do i=1,nc-1
c        p    = gam1*( q(3,i) - 0.5d0*q(2,i)**2/q(1,i) )
         p1   = gam1*( q(3,i) - 0.5d0*q(2,i)**2/q(1,i) )
         p2   = gam1*( q(3,i+1) - 0.5d0*q(2,i+1)**2/q(1,i+1) )

c        Integral of pressure
c        cost = cost + p * dx
         cost = cost + 0.5d0*(p1+p2) * dx

c        Pressure at some cell: WARNING if you change grid you must change
c        the cell number
c        if(i.eq.28) cost = p
      enddo

      return
      end

C Cost function, L2 norm of pressure difference
C Uses P1 interpolation of cell-center pressures
C Extra work to calculate integral in first and last cell
      subroutine costfunc2(nc, q, cost)
      implicit none
      include 'param.inc'
      integer          nc
      double precision q(3,*), cost

      double precision p
      integer          i

      cost = cost + pinf*dx/4.0d0 + pout*dx/4.0d0

      i    = 1
      p    = gam1*( q(3,i) - 0.5d0*q(2,i)**2/q(1,i) )
      cost = cost + p*(dx/4.0d0 + dx/2.0)

      do i=2,nc-1
         p    = gam1*( q(3,i) - 0.5d0*q(2,i)**2/q(1,i) )
         cost = cost + p * dx
      enddo

      i    = nc
      p    = gam1*( q(3,i) - 0.5d0*q(2,i)**2/q(1,i) )
      cost = cost + p*(dx/4.0d0 + dx/2.0)

      return
      end

C Returns nozzle cross-section area
      double precision function nozzle(x)
      implicit none
      include 'param.inc'
      double precision x
      ! S-shaped diverging nozzle
      if(noztyp.eq.1)then
         nozzle = da + db*dtanh(dc*x - dd)
      ! Converging-diverging nozzle: This is not C^2
      elseif(noztyp.eq.2)then
         if(x.lt.-0.5d0 .or. x.gt.0.5d0)then
            nozzle = 2.0d0
         else
            nozzle = 1.0d0 + dsin(M_PI*x)**2
         endif
      ! Straight diverging nozzle
      elseif(noztyp.eq.3)then
         nozzle = ain + x*(aout-ain)/L
      ! Polynomial converging-diverging nozzle
      elseif(noztyp.eq.4)then
         if(x.lt.-1.0d0 .or. x.gt.+1.0d0)then
            nozzle = ain
         else
            nozzle = ath + (ain-ath)*(3.0d0*x**2 - 3.0d0*x**4 + x**6)
         endif
      endif
      return
      end

C Maximum of two numbers
      double precision function dmax1(a, b)
      implicit none
      double precision a, b

      if(a .gt. b) then
         dmax1 = a
      else
         dmax1 = b
      endif

      return
      end

C Minimum of two numbers
      double precision function dmin1(a, b)
      implicit none
      double precision a, b

      if(a .lt. b) then
         dmin1 = a
      else
         dmin1 = b
      endif

      return
      end

C.....Error function, from Abromovitz and Stegun
      double precision function ERRF(X)
      double precision X,ARG,E,VB,T,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6

      ARG = X*X
      if(ARG .lt. 20.0d0)then
            E = dexp(-ARG)
      else
            E = 0.0d0
      endif
      VB = dabs(X)
      T = 1.0d0/(1.0d0 + 0.3275911d0*VB)
      tmp1 = 1.061405429d0*T
      tmp2 = (tmp1 - 1.453152027d0)*T
      tmp3 = (tmp2 + 1.421413741d0)*T
      tmp4 = (tmp3 - 0.284496736d0)*T
      tmp5 = (tmp4 + 0.254829592d0)*T
      tmp6 = 1.0d0 - tmp5*E
      if(X .lt. 0.0d0)then
            ERRF = -tmp6
      else
            ERRF =  tmp6
      endif
      return
      end
c
c     MUSCL reconstruction
c
      subroutine muscl(qjm1, qj, qjp1, qjph)
      implicit none
      include 'param.inc'
      double precision qjm1(*), qj(*), qjp1(*), qjph(*)

      integer i
      double precision dl, dr, epsq, d1, d2, xm, kkk

      epsq = 1.0d-8
      kkk  = 1.0d0/3.0d0

      if(tolimit.eq.0)then

         do i=1,3
            qjph(i) = qj(i)
         enddo

      else

      do i=1,3
         d1 = qj(i)   - qjm1(i)
         d2 = qjp1(i) - qj(i)

         XM = ( D1*(D2*D2+2.E0*EPSQ)
     1        + D2*(2.E0*D1*D1+EPSQ) )
     2                / ( 2.E0*D1*D1 - D1*D2 + 2.E0*D2*D2 + 3.E0*EPSQ )

         qjph(i) = qj(i) + 0.5d0 * xm
c        Unlimited case
c        qjph(i) = qj(i) + 0.25d0 * (d1 + d2)
      enddo

      endif

      return
      end
c
c     Roe flux function
c
      subroutine Flux_Roe(a, ql, qr, resl, resr)
      implicit none
      include 'param.inc'
      double precision a, ql(3), qr(3), resl(3), resr(3), f(3)
      double precision rl, ul, pl, hl, cl
      double precision rr, ur, pr, hr, cr
      double precision ra, ua, pa, ha, ca
      double precision sqrl, sqrr
      double precision dq(3), m1, m2, a1, a2, a3
      double precision l1, l2, l3
      double precision EPS

      EPS = 1.0d-2

      rl = ql(1)
      ul = ql(2)/ql(1)
      pl = (gam-1.0d0)*(ql(3) - 0.5d0*rl*ul**2)
      cl = dsqrt(gam*pl/rl)
      hl = cl**2/(gam-1.0d0) + 0.5d0*ul**2

      rr = qr(1)
      ur = qr(2)/qr(1)
      pr = (gam-1.0d0)*(qr(3) - 0.5d0*rr*ur**2)
      cr = dsqrt(gam*pr/rr)
      hr = cr**2/(gam-1.0d0) + 0.5d0*ur**2

      sqrl = dsqrt(rl)
      sqrr = dsqrt(rr)
      ua = (sqrl*ul + sqrr*ur)/(sqrl + sqrr)
      ha = (sqrl*hl + sqrr*hr)/(sqrl + sqrr)
      ca = dsqrt((gam-1.0d0)*(ha - 0.5d0*ua**2))

      dq = qr - ql

      m1 = (dq(2) - ua*dq(1))/ca
      m2 = (gam-1.0d0)*(dq(3) + 0.5d0*ua**2*dq(1) - ua*dq(2))/ca**2
      a3 = 0.5d0*(m1 + m2)
      a1 = 0.5d0*(m2 - m1)
      a2 = dq(1) - a1 - a3

      l1 = dabs(ua-ca)
      l2 = dabs(ua)
      l3 = dabs(ua+ca)

c     entropy fix
c     if(l1.lt.EPS) l1 = 0.5d0*(EPS + l1*l1/EPS)

c     central flux
      f(1) = 0.5d0*( rl*ul + rr*ur )
      f(2) = 0.5d0*( pl+rl*ul**2 + pr + rr*ur**2 )
      f(3) = 0.5d0*( rl*hl*ul + rr*hr*ur )

      f(1) = f(1) - 0.5d0*( a1*l1 + a2*l2 + a3*l3 )
      f(2) = f(2) - 0.5d0*( a1*l1*(ua-ca) + a2*l2*ua + a3*l3*(ua+ca) )
      f(3) = f(3) - 0.5d0*( a1*l1*(ha-ua*ca) + a2*l2*ua**2/2.0d0 +
     1                    a3*l3*(ha+ua*ca) )

      resl(:) = resl(:) + a*f(:)
      resr(:) = resr(:) - a*f(:)

      return
      end
