      integer rk4, rktvd
      parameter(rk4=1, rktvd=2)

      integer ictype, iorder, inttyp, rkmeth
      real kkk
      parameter(kkk=1.0/3.0)
      common/iparams/ictype,iorder,inttyp,rkmeth

      real ark(3), brk(3)
      common/rkparams/ark, brk

      integer LINEAR, BURGER
      parameter(LINEAR=1, BURGER=2)
      integer ieq
      common/equation/ieq

      real wpower
      common/weight/wpower
