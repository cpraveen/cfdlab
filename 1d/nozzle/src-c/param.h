#ifdef DEFINE_GLOBALS
#define GLOBAL
#else
#define GLOBAL extern
#endif

      GLOBAL double gam, gam1, gasconst, cfl;

      GLOBAL double dx, ain, aout, da, db, dc, dd, L;

      GLOBAL double rinf, uinf, pinf, minf, Tinf, pout, prat;

      GLOBAL int maxiter, flux1, flux2;

