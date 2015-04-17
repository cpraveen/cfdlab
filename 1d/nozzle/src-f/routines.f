

C Flux for inflow face
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
      subroutine costfunc(q, ptarg, cost)
      implicit none
      include 'param.inc'
      double precision q(3), ptarg, cost

      double precision p

      p    = gam1*( q(3) - 0.5d0*q(2)**2/q(1) )
      cost = cost + ( p - ptarg )**2

      return
      end

C Returns nozzle cross-section area
      double precision function nozzle(L, ain, aout, x, c, d)
      implicit none
      double precision L, ain, aout, x, a, b, c, d
      b  = (aout - ain)/(dtanh(L*c - d) - dtanh(-d))
      a  = ain - b*dtanh(-d)
      nozzle = a + b*dtanh(c*x - d)
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
