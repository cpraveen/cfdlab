      integer npmax, ntmax, nnbrmax, nxmax, nymax
      parameter(nxmax=201, nymax=201)
      parameter(npmax=nxmax*nymax, ntmax=2*npmax, nnbrmax=20)
      integer CONSTANT, INTERIOR, PERIODIC
      parameter(CONSTANT=0, INTERIOR=1, PERIODIC=2)

      double precision coord(2,npmax), var(npmax), var_old(npmax), 
     +                 ax(nnbrmax,npmax), ay(nnbrmax,npmax), 
     +                 res(npmax), varx(npmax), vary(npmax)
      double precision xmin, xmax, ymin, ymax, dx, dy, lx, ly
      double precision dt, Ttime, kkk
      double precision varmin0, varmax0
      double precision pertx, perty
      common/dvars/coord,var,var_old,ax,ay,res,varx,vary,xmin,xmax,ymin,
     +             ymax,dx,dy,lx,ly,dt,Ttime,kkk,varmin0,varmax0,
     +             pertx,perty

      integer          nnum(0:nxmax,0:nymax), nnbr(npmax), ptype(npmax),
     +                 conn(nnbrmax,npmax), nx, ny, tri(3,ntmax)
      integer          np0, npts, ntri, veltype, ictype, order
      integer          makeanim, animinterval, perturb
      common/ivars/nnum,nnbr,ptype,conn,nx,ny,tri,np0,npts,ntri,veltype,
     +             ictype,order,makeanim,animinterval,perturb
