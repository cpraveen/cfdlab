      integer ifirst, itvd, iweno3, iweno5, recon_scheme
      parameter(ifirst=0,itvd=1,iweno3=2,iweno5=3)
      common/scheme/recon_scheme

      real w1
      common/poslim/w1

      real ark(3), brk(3)
      common/ssprk/ark,brk
