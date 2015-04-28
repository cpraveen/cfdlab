      integer ifirst, itvd, iweno, recon_scheme
      parameter(ifirst=0,itvd=1,iweno=2)
      common/scheme/recon_scheme

      real    mu
      parameter(mu=2.0)

      real ark(3), brk(3)
      common/ssprk/ark,brk
