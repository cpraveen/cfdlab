      integer nel
      parameter(nel=100)
      real u(2,nel), v(2,nel), uold(2,nel), vold(2,nel)
      real ua(nel), uaold(nel)
      real ures(2,nel), vres(2,nel), ux(nel), uxold(nel)
      real B(nel+1)
      real x(2,nel)
      real cfl
      common/sol/u,v,uold,vold,ures,vres,ux,uxold,B,x,cfl
