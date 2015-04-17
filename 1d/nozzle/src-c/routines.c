#include <math.h>
#include "double.h"
#include "param.h"

#include <stdio.h>
#include <stdlib.h>
#include "routines.h"

/* Computes finite volume residual */
void residu(int nc, int nf, REALX * a, REALQ ** q, REALA ** res)
{
   int i, j, nl, nr;

   for(i = 0; i < nc; i++)
      for(j = 0; j < 3; j++)
         res[i][j] = 0.0;

   flux_in(a[0], res[0]);
   for(i = 1; i < nf - 1; i++) {
      nl = i - 1;
      nr = i;
      if(flux1 == 1)
         flux_ausm(a[i], q[nl], q[nr], res[nl], res[nr]);
#ifndef NO_KFVS_
      else if(flux1 == 2)
         flux_kfvs(a[i], q[nl], q[nr], res[nl], res[nr]);
#endif
      else if(flux1 == 3)
         flux_lf(a[i], q[nl], q[nr], res[nl], res[nr]);
      else {
         printf("Flux type is not known, %d\n", flux1);
         exit(0);
      }
   }
   flux_out(a[nf - 1], q[nf - 2], q[nf - 3], res[nf - 2]);

   for(i = 0; i < nc; i++)
      source_term(a[i], a[i + 1], q[i], res[i]);

}

//Flux for inflow face
void flux_in(REALX a, REALA * res)
{
   int i;
   REALA einf, fl[3];

   fl[0] = rinf * uinf;
   fl[1] = pinf + rinf * uinf * uinf;
   einf = pinf / gam1 + 0.5 * rinf * uinf * uinf;
   fl[2] = (einf + pinf) * uinf;

   for(i = 0; i < 3; i++)
      res[i] = res[i] - a * fl[i];
}

//Flux for interior face
//AUSMDV flux, taken from nozzle code of Manoj Nair, CTFD, NAL
void flux_ausm(REALX a, REALQ * ql, REALQ * qr, REALA * resl, REALA * resr)
{
   int i;
   REALA rho_left, rho_right, u_left, u_right,
      p_left, p_right, c_left, c_right, c_mean,
      H_left, H_right, amach_l, amach_r,
      alpha_l, alpha_r, temp1, temp2, temp3, u_lp, u_rm,
      p_lp, p_rm, ak_half, p_half, s_half, flux_v1,
      flux_d1, flux_p, flux_n, fl[3], p_min;

   rho_left = ql[0];
   u_left = ql[1] / ql[0];
   p_left = gam1 * (ql[2] - 0.5 * ql[1] * ql[1] / ql[0]);

   rho_right = qr[0];
   u_right = qr[1] / qr[0];
   p_right = gam1 * (qr[2] - 0.5 * qr[1] * qr[1] / qr[0]);

   c_left = sqrt(gam * p_left / rho_left);
   H_left = c_left * c_left / gam1 + 0.5 * u_left * u_left;

//......Right State 
   c_right = sqrt(gam * p_right / rho_right);
   H_right = c_right * c_right / gam1 + 0.5 * u_right * u_right;

//......Mean Speed of Sound
   c_mean = (c_left > c_right) ? c_left : c_right;

//......Left and Right Mach no.
   amach_l = fabs(u_left) / c_mean;
   amach_r = fabs(u_right) / c_mean;

//......ALPHA
   temp1 = p_left / rho_left;
   temp2 = p_right / rho_right;

   alpha_l = 2 * temp1 / (temp1 + temp2);
   alpha_r = 2 * temp2 / (temp1 + temp2);

//......U_LEFT+
   temp1 = u_left + c_mean;
   temp2 = u_left + fabs(u_left);
   if(amach_l <= 1)
      u_lp = alpha_l * (0.25 * temp1 * temp1 / c_mean - 0.5 * temp2)
         + 0.5 * temp2;
   else
      u_lp = 0.5 * temp2;

//......U_RIGHT-
   temp1 = u_right - c_mean;
   temp2 = u_right - fabs(u_right);
   if(amach_r <= 1)
      u_rm = alpha_r * (-0.25 * temp1 * temp1 / c_mean - 0.5 * temp2)
         + 0.5 * temp2;
   else
      u_rm = 0.5 * temp2;

//......P_LEFT+
   temp1 = u_left / c_mean + 1;
   temp2 = 2 - u_left / c_mean;
   if(amach_l <= 1)
      p_lp = 0.25 * p_left * temp1 * temp1 * temp2;
   else
      p_lp = 0.5 * p_left * (u_left + fabs(u_left)) / (u_left);

//......P_RIGHT-
   temp1 = u_right / c_mean - 1;
   temp2 = 2 + u_right / c_mean;
   if(amach_r <= 1)
      p_rm = 0.25 * p_right * temp1 * temp1 * temp2;
   else
      p_rm = 0.5 * p_right * (u_right - fabs(u_right)) / (u_right);

   ak_half = 10.0;
   p_half = p_lp + p_rm;
   p_min = (p_left < p_right) ? p_left : p_right;
   s_half = fabs(p_right - p_left) / p_min;
   temp3 = ak_half * s_half;
   if(temp3 < 1.0)
      s_half = 0.5 * temp3;
   else
      s_half = 0.5;

//......Calculating the inviscid interface fluxes

//......Mass flux
   fl[0] = u_lp * rho_left + u_rm * rho_right;

//......Enthalpy flux
   fl[2] = 0.5 * (fl[0] * (H_left + H_right)
                  - fabs(fl[0]) * (H_right - H_left));

//......Momentum flux

   flux_v1 = rho_left * u_left * u_lp + rho_right * u_right * u_rm;

   flux_d1 = 0.5 * (fl[0] * (u_left + u_right) -
                    fabs(fl[0]) * (u_right - u_left));

   flux_p = p_half;

   flux_n = (0.5 + s_half) * flux_v1 + (0.5 - s_half) * flux_d1;

   fl[1] = flux_n + flux_p;

   for(i = 0; i < 3; i++) {
      resl[i] = resl[i] + a * fl[i];
      resr[i] = resr[i] - a * fl[i];
   }

}

//Flux for interior face
//Lax-Fredrichs flux
void flux_lf(REALX a, REALQ * ql, REALQ * qr, REALA * resl, REALA * resr)
{
   int i;
   REALA rho_left, rho_right, u_left, u_right,
      p_left, p_right, c_left, c_right,
      h_left, h_right, rl, rr, rs, ua, ha, aa, lam, fl[3], fr[3], flux[3];

   rho_left = ql[0];
   u_left = ql[1] / ql[0];
   p_left = gam1 * (ql[2] - 0.5 * ql[1] * ql[1] / ql[0]);
   c_left = sqrt(gam * p_left / rho_left);
   h_left = c_left * c_left / gam1 + 0.5 * u_left * u_left;


//......Right State 
   rho_right = qr[0];
   u_right = qr[1] / qr[0];
   p_right = gam1 * (qr[2] - 0.5 * qr[1] * qr[1] / qr[0]);
   c_right = sqrt(gam * p_right / rho_right);
   h_right = c_right * c_right / gam1 + 0.5 * u_right * u_right;

   fl[0] = ql[1];
   fl[1] = p_left + rho_left * u_left * u_left;
   fl[2] = ql[1] * h_left;

   fr[0] = qr[1];
   fr[1] = p_right + rho_right * u_right * u_right;
   fr[2] = qr[1] * h_right;

//......Roe averages
   rl = sqrt(rho_left);
   rr = sqrt(rho_right);
   rs = 1.0 / (rl + rr);
   ua = (rl * u_left + rr * u_right) * rs;
   ha = (rl * h_left + rr * h_right) * rs;
   aa = sqrt((gam - 1.0) * (ha - 0.5 * ua * ua));

//    ua = 0.5*( u_left + u_right);
//    aa = 0.5*( c_left + c_right);

   lam = fabs(ua) + aa;

   for(i = 0; i < 3; i++) {
      flux[i] = 0.5 * (fl[i] + fr[i]) - 0.5 * lam * (qr[i] - ql[i]);
      resl[i] = resl[i] + a * flux[i];
      resr[i] = resr[i] - a * flux[i];
   }
}

#ifndef NO_KFVS_
//Flux for interior face
//KFVS flux, taken from nozzle code of Manoj Nair
void flux_kfvs(REALX a, REALQ *ql, REALQ *qr, REALA *resl, REALA *resr)
{
   int i;
   REALA rl, rr, ul, ur, pl, pr, el, er, betal, betar,
      sl, sr, Al, Ar, Bl, Br, Fp[3], Fm[3];

   rl = ql[0];
   ul = ql[1] / ql[0];
   pl = gam1 * (ql[2] - 0.5 * ql[1] * ql[1] / ql[0]);
   el = pl / gam1 + 0.5 * rl * ul * ul;

   rr = qr[0];
   ur = qr[1] / qr[0];
   pr = gam1 * (qr[2] - 0.5 * qr[1] * qr[1] / qr[0]);
   er = pr / gam1 + 0.5 * rr * ur * ur;

   betal = 0.5 * rl / pl;
   sl = ul * sqrt(betal);
   Al = 0.5 * (1.0 + erf(sl));
   Bl = 0.5 * exp(-sl * sl) / sqrt(M_PI * betal);

   Fp[0] = rl * (ul * Al + Bl);
   Fp[1] = (pl + rl * ul * ul) * Al + rl * ul * Bl;
   Fp[2] = (el + pl) * ul * Al + (el + 0.5 * pl) * Bl;

//    Negative flux
   betar = 0.5 * rr / pr;
   sr = ur * sqrt(betar);
   Ar = 0.5 * (1.0 - erf(sr));
   Br = 0.5 * exp(-sr * sr) / sqrt(M_PI * betar);

   Fm[0] = rr * (ur * Ar - Br);
   Fm[1] = (pr + rr * ur * ur) * Ar - rr * ur * Br;
   Fm[2] = (er + pr) * ur * Ar - (er + 0.5 * pr) * Br;

   for(i = 0; i < 3; i++) {
      resl[i] = resl[i] + a * (Fp[i] + Fm[i]);
      resr[i] = resr[i] - a * (Fp[i] + Fm[i]);
   }
}
#endif

//Flux for outflow face
void flux_out(REALX a, REALQ * ql, REALQ * qll, REALA * res)
{
   int i;
   REALA rl, ul, rll, ull;
   REALA r, u, e, fl[3];

   rl = ql[0];
   ul = ql[1] / ql[0];

   rll = qll[0];
   ull = qll[1] / qll[0];

   r = 1.5 * rl - 0.5 * rll;
   u = 1.5 * ul - 0.5 * ull;
   e = pout / gam1 + 0.5 * r * u * u;

   fl[0] = r * u;
   fl[1] = pout + r * u * u;
   fl[2] = (e + pout) * u;

   for(i = 0; i < 3; i++)
      res[i] = res[i] + a * fl[i];
}

//Source term contribution to residual
void source_term(REALX al, REALX ar, REALQ * q, REALA * res)
{
   REALA p;

   p = gam1 * (q[2] - 0.5 * q[1] * q[1] / q[0]);
   res[1] = res[1] - p * (ar - al);
}

//Cost function, L2 norm of pressure difference
void costfunc(REALQ * q, double ptarg, REALA * cost)
{
   REALA p, dp;

   p = gam1 * (q[2] - 0.5 * q[1] * q[1] / q[0]);
   dp = p - ptarg;
   *cost = *cost + dp * dp;
}
