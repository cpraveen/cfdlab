      double precision gam, gam1, gasconst, cfl
      common/gas/gam, gam1, gasconst, cfl

      double precision dx, ain, aout, da, db, dc, dd, L
      common/grid/dx, ain, aout, da, db, dc, dd, L

      double precision rinf, uinf, pinf, minf, Tinf, pout, prat
      common/flow/rinf, uinf, pinf, minf, Tinf, pout, prat

      double precision M_PI
      parameter(M_PI=3.14159265358979323846d0)

      integer maxiter, flux1, flux2
      common/param/maxiter,flux1,flux2
