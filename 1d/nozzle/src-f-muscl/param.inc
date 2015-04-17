      double precision gam, gam1, gasconst, cfl
      common/gas/gam, gam1, gasconst, cfl

      double precision dx, ain, ath, aout, da, db, dc, dd, L, xmin
      common/grid/dx, ain, ath, aout, da, db, dc, dd, L, xmin

      double precision rinf, uinf, pinf, minf, Tinf, pout, prat, qinf(3)
      common/flow/rinf, uinf, pinf, minf, Tinf, pout, prat, qinf

      double precision ark(3)
      common/time/ark

      double precision M_PI
      parameter(M_PI=3.14159265358979323846d0)

      double precision kappa
      parameter(kappa=0.3333333333333333d0)

      integer maxiter, flux1, flux2, noztyp, tolimit
      common/param/maxiter,flux1,flux2,noztyp,tolimit
