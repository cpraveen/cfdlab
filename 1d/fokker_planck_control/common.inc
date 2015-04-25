      integer ifirst, itvd, iweno, recon_scheme
      parameter(ifirst=0,itvd=1,iweno=2)
      common/scheme/recon_scheme

      real    mu
      parameter(mu=2.0)

      integer nc
      parameter(nc=500)

      real xmin, xmax
      parameter(xmin=-2.0, xmax=2.0)

      real nu_reg
      parameter(nu_reg=0.25)

      integer initc
      common/init/initc
