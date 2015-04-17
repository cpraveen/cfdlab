c------------------------------------------------------------------------------
c     Solves linear advection equation on 1-d overset grids
c------------------------------------------------------------------------------
      program solve
      implicit none
      include 'param.h'
      include 'vars.h'
      include 'common.h'

      integer ifid, i, j
      integer niter, nit, irk, nrk
      real    cfl, dt, dtstage, time, dxmin
      real    errm, err2

c     Equation: linear advection or burgers equation
      ieq = LINEAR

c     Initial condition
c     1=gaussian, 2=step
      ictype = 1

c     Order: 1=first order, 2=second order, 3=limited second order
      iorder = 2

c     Interpolation meth: 1=linear, 2=quadratic, 3=cubic, 4=meshless
c     5=3'rd order meshless
      inttyp = 4

c     weight power for meshless method
      wpower = 1.0
      
c     CFL number
      cfl   = 0.95

c     Number of time steps
      niter = 111

c     RK time integration scheme: rk4, rktvd
c     rkmeth = rk4
      rkmeth = rktvd

c     Number of rk stages
      if(rkmeth.eq.rk4) nrk = 4
      if(rkmeth.eq.rktvd) nrk = 3

c     Coefficients for TVD RK scheme
      ark(1) = 0.0
      ark(2) = 3.0/4.0
      ark(3) = 1.0/3.0

      brk(1) = 1.0 - ark(1)
      brk(2) = 1.0 - ark(2)
      brk(3) = 1.0 - ark(3)

c     Read grids from file
      print*,'Reading grids'
      ifid=10
      open(ifid, file='grids.dat')
      read(ifid,*) ngrid
      read(ifid,*)(np(i), i=1,ngrid)
      read(ifid,*)(dx(i), i=1,ngrid)
      read(ifid,*)((x(i,j),utype(i,j),j=1,np(i)), i=1,ngrid)
      close(ifid)

      print*,'Number of grids =',ngrid
      write(*,*)(np(i), i=1,ngrid)

c     Read overlap information
      print*,'Reading overlap info from overlap.dat'
      open(ifid, file='overlap.dat')
      read(ifid,*) nopts
      do i=1,nopts
         read(ifid,*)(olap(j,i),j=1,7)
      enddo
      close(ifid)

c     Enlarge stencil if required, but only for second order meshless
      if(inttyp.eq.4)then
         call stable_stencil()
      endif

#ifdef _GRAFIC
      CALL GRINIT(5,6,'DEMONSTRATION OF GRAFIC')
      CALL GRSSET(0.0, 2.0, 0.0, 1.0)
#endif

c     Set initial condition
      call SetInitCond

#ifdef _GRAFIC
      call PlotGrafic(0.0)
#endif

c     Save initial condition
      call SaveSol
      call system("mv -f sol.dat sol0.dat")

c     Time step
      dxmin = 1.0e20
      do i=1,ngrid
         dxmin = min(dxmin, dx(i))
      enddo
      dt = cfl*dxmin
      
      niter = 1.0/dt

      print*,'Time step =',dt

c     Loop over time iterations
      print*,'Starting time integration'
      time = 0.0
      open(55, file='hist.dat')
      do nit=1,niter

c        Store solution at previous time level
         var_old(:,:) = var(:,:)

c        Loop over RK stages
         do irk=1,nrk

            if(rkmeth.eq.rk4)then
               dtstage = dt/(nrk-irk+1)
            else
               dtstage = dt
            endif

c           Store intermediate stage values
c           USE var0 FOR RHS COMPUTATION
c           USE var FOR INTERPOLATION
            var0(:,:) = var(:,:)

c           Update grid points
            do i=1,ngrid
               call UpdateGrid(irk, i, dtstage)
            enddo

c           Update overlap points
            call UpdateOverlap(irk, dtstage)

         enddo

         time = time + dt

#ifdef _GRAFIC
         call PlotGrafic(time)
#endif

         call Errors(time, errm, err2)
         print*,nit,time,errm,err2
         write(55,*)nit,time,errm,err2

      enddo
      close(55)

c     Solve the solution
      call SaveSol

      stop
      end
c------------------------------------------------------------------------------
c     Enlarge stencil for better stability
c------------------------------------------------------------------------------
      subroutine stable_stencil()
      include 'param.h'
      include 'vars.h'
      integer i, irg, irp, idg, id1, id2
      real    ds, fact

      print*,'Setting stable stencil'

      do i=1,nopts
c        Receiving grid and point
         irg = olap(1,i)
         irp = olap(2,i)

c        Donor grid and point
         idg = olap(3,i)
         id1 = olap(4,i)
         id2 = olap(5,i)

         ds  =   x(idg,id2) - x(idg,id1)
         fact= ( x(irg,irp) - x(idg,id1) )/ds

         if(fact.le.0.5)then
            olap(4,i) = id1 - 1
         else if(fact.gt.0.5)then
            olap(5,i) = id2 + 1
         endif
      enddo

      end
c------------------------------------------------------------------------------
c     Set initial condition
c------------------------------------------------------------------------------
      subroutine SetInitCond
      implicit none
      include 'param.h'
      include 'vars.h'

      integer i, j
      real    InitCond

      print*,'Setting initial condition'
      do i=1,ngrid
         do j=1,np(i)
            var(i,j) = InitCond(x(i,j))
         enddo
      enddo

      var0(:,:) = var(:,:)
      return
      end
c------------------------------------------------------------------------------
c     Function returns initial condition
c------------------------------------------------------------------------------
      real function InitCond(xx)
      implicit none
      include 'common.h'
      real xx

      real dx

      if(ictype.eq.1)then
         dx       = xx - 0.5
         InitCond = exp(-100.0*dx**2)
      else if(ictype.eq.2)then
         if(xx.lt.0.5)then
            InitCond = 1.0
         else
            InitCond = 0.0
         endif
      endif

      return
      end
c------------------------------------------------------------------------------
c     Update solution on igrid
c------------------------------------------------------------------------------
      subroutine UpdateGrid(irk, igrid, dt)
      implicit none
      include 'param.h'
      include 'vars.h'
      include 'common.h'
      integer irk, igrid
      real    dt
      real    fluxfunc

      integer i
      real    ul, ur, fl, fr, muscl

      do i=1,np(igrid)
         if(utype(igrid,i).eq.UPDATE)then
c           Left face flux
            ul = var0(igrid,i-1) +
     1           muscl( var0(igrid,i-2),
     1                  var0(igrid,i-1),
     2                  var0(igrid,i  ) )
            ur = var0(igrid,i) -
     1           muscl( var0(igrid,i-1),
     1                  var0(igrid,i),
     2                  var0(igrid,i+1) )
            fl = fluxfunc(ul,ur)
c           Right face flux
            ul = var0(igrid,i) +
     1           muscl( var0(igrid,i-1),
     1                  var0(igrid,i),
     2                  var0(igrid,i+1) )
            ur = var0(igrid,i+1) -
     1           muscl( var0(igrid,i),
     1                  var0(igrid,i+1),
     2                  var0(igrid,i+2) )
            fr = fluxfunc(ul,ur)
c           Perform RK update
            if(rkmeth.eq.rk4)then
               var(igrid,i) = var_old(igrid,i) - dt*(fr - fl)/dx(igrid)
            elseif(rkmeth.eq.rktvd)then
               var(igrid,i) = ark(irk)*var_old(igrid,i) +
     1                        brk(irk)*(var0(igrid,i) - 
     2                          dt*(fr - fl)/dx(igrid) )
            else
               print*,'Uknown RK scheme'
               stop
            endif
         endif
      enddo

      return
      end
c------------------------------------------------------------------------------
c     MUSCL reconstruction
c------------------------------------------------------------------------------
      real function muscl(ul, u, ur)
      implicit none
      include 'common.h'
      real ul, u, ur

      real limiter, koren_delta, albada_delta
      real dl, dr, s

      dl  = u - ul
      dr  = ur - u
      s   = limiter(dl, dr)    

c     muscl = u + 0.25*s*( (1.0-kkk)*dl + (1.0+kkk)*dr )
c     muscl = u + 0.5d0*koren_delta(dl, dr)
c     muscl = u + 0.5d0*albada_delta(dl, dr)

      muscl = 0.25*s*( (1.0-kkk)*dl + (1.0+kkk)*dr )
c     muscl = 0.5*koren_delta(dl, dr)

      return
      end
c------------------------------------------------------------------------------
c------------------------------------------------------------------------------
      real function limiter(a, b)
      implicit none
      include 'common.h'
      real a, b

      real eps, num, den

      if(iorder.eq.1)then
         limiter = 0.0
         return
      else if(iorder.eq.2)then
         limiter = 1.0
         return
      else if(iorder.eq.3)then
         eps = 1.0e-8
         num = 3.0*a*b + eps
         den = 2.0*(a-b)**2 + num
         limiter = num/den
         return
      else if(iorder.eq.4)then
         eps = 1.0e-8
         num = 2.0*a*b + eps
         den = a**2 + b**2 + eps
         limiter = num/den
         return
      else
         print*,'limiter: unknown iorder =',iorder
         stop
      endif

      return
      end
c------------------------------------------------------------------------------
c     Delta form of koren limiter
c------------------------------------------------------------------------------
      real function koren_delta(a, b)
      implicit none
      include 'common.h'
      real a, b

      real eps, num, den

      eps = 1.0e-8
      num = (2.0*a*a + eps)*b + (b*b + 2.0*eps)*a
      den = 2.0*a*a + 2.0*b*b - a*b + 3.0*eps
      koren_delta = num/den

      return
      end
c------------------------------------------------------------------------------
c     Delta form of van albada limiter
c------------------------------------------------------------------------------
      real function albada_delta(a, b)
      implicit none
      include 'common.h'
      real a, b

      real eps, num, den

      eps = 1.0e-8
      num = (a*a + eps)*b + (b*b + eps)*a
      den = a*a + b*b + 2.0*eps
      albada_delta = num/den

      return
      end
c------------------------------------------------------------------------------
c     Numerical flux function
c------------------------------------------------------------------------------
      real function fluxfunc(ul, ur)
      implicit none
      include 'common.h'
      real ul, ur

      real a

      if(ieq.eq.LINEAR)then
         fluxfunc = ul
      else if(ieq.eq.BURGER)then
         if(ul.eq.ur)then
            fluxfunc = 0.5*ul*ul
         else
            a = 0.5*(ul + ur)
            if(a.gt.0.0)then
               fluxfunc = 0.5*ul*ul
            else
               fluxfunc = 0.5*ur*ur
            endif
         endif
      else
         print*,'Unknown equation'
         stop
      endif

      return
      end
c------------------------------------------------------------------------------
c     Update points in overlap region
c------------------------------------------------------------------------------
      subroutine UpdateOverlap(irk, dt)
      implicit none
      include 'common.h'
      integer irk
      real    dt

      if(inttyp.eq.1)then
         call UpdateOverlap_Interp_Linear
      else if(inttyp.eq.2)then
         call UpdateOverlap_Interp_Quadratic
      else if(inttyp.eq.3)then
         call UpdateOverlap_Interp_Cubic
      else if(inttyp.eq.4)then
         call UpdateOverlap_meshless(irk, dt)
      else if(inttyp.eq.5)then
         call UpdateOverlap_meshless3(irk, dt)
      else
         print*,'UpdateOverlap: unknown inttyp =',inttyp
         stop
      endif

      return
      end
c------------------------------------------------------------------------------
c     Update points in overlap region using linear interpolation
c------------------------------------------------------------------------------
      subroutine UpdateOverlap_Interp_Linear
      implicit none
      include 'param.h'
      include 'vars.h'

      integer i, irg, irp, idg, id1, id2
      real    ds, fact

      do i=1,nopts
c        Receiving grid and point
         irg = olap(1,i)
         irp = olap(2,i)

c        Donor grid and point
         idg = olap(3,i)
         id1 = olap(4,i)
         id2 = olap(5,i)

         ds  =   x(idg,id2) - x(idg,id1)
         fact= ( x(irg,irp) - x(idg,id1) )/ds

         var(irg,irp) = (1.0-fact)*var(idg,id1) + fact*var(idg,id2)
      enddo

      return
      end
c------------------------------------------------------------------------------
c     Update points in overlap region quadratic interpolation
c------------------------------------------------------------------------------
      subroutine UpdateOverlap_Interp_Quadratic
      implicit none
      include 'param.h'
      include 'vars.h'

      integer i, j, k, irg, irp, idg, id(4), nn
      real    xp, vp, xx(4), Lagrange

      nn = 3

      do i=1,nopts
c        Receiving grid and point
         irg = olap(1,i)
         irp = olap(2,i)

c        Donor grid and point
         idg = olap(3,i)
         id(1) = olap(4,i)
         id(2) = olap(5,i)
         id(3) = olap(6,i)
         id(4) = olap(7,i)

         xp = x(irg,irp)
         do j=1,nn
            xx(j) = x(idg, id(j))
         enddo

         vp = 0.0
         do j=1,nn
            vp = vp + Lagrange(nn, j, xx, xp)*var(idg, id(j))
         enddo

         var(irg,irp) = vp

      enddo

      return
      end
c------------------------------------------------------------------------------
c     Update points in overlap region using cubic interpolation
c------------------------------------------------------------------------------
      subroutine UpdateOverlap_Interp_Cubic
      implicit none
      include 'param.h'
      include 'vars.h'

      integer i, j, k, irg, irp, idg, id(4), nn
      real    xp, vp, xx(4), Lagrange

      nn = 4

      do i=1,nopts
c        Receiving grid and point
         irg = olap(1,i)
         irp = olap(2,i)

c        Donor grid and point
         idg = olap(3,i)
         id(1) = olap(4,i)
         id(2) = olap(5,i)
         id(3) = olap(6,i)
         id(4) = olap(7,i)

         xp = x(irg,irp)
         do j=1,nn
            xx(j) = x(idg, id(j))
         enddo

         vp = 0.0
         do j=1,nn
            vp = vp + Lagrange(4, j, xx, xp)*var(idg, id(j))
         enddo

         var(irg,irp) = vp

      enddo

      return
      end
c------------------------------------------------------------------------------
c------------------------------------------------------------------------------
      real function Lagrange(n, i, xx, xp)
      implicit none
      integer n, i
      real    xx(*), xp

      integer j

      Lagrange = 1.0
      do j=1,n
         if(j.ne.i)then
            Lagrange = Lagrange * ( xp - xx(j) ) / ( xx(i) - xx(j) )
         endif
      enddo

      return
      end
c------------------------------------------------------------------------------
c     Update points in overlap region using meshless method
c------------------------------------------------------------------------------
      subroutine UpdateOverlap_meshless(irk, dt)
      implicit none
      include 'param.h'
      include 'vars.h'
      include 'common.h'
      integer irk
      real    dt

      integer i, irg, irp, idg, id1, id2
      real    vx1, vx2, vxp, dx1, dx2, dnm, f1, f2, fp, res, dl, dr
      real    ul, ur, u1, u2, s, limiter, w1, w2
      real    fluxfunc

      do i=1,nopts
c        Receiving grid and point
         irg = olap(1,i)
         irp = olap(2,i)

c        Donor grid and point
         idg = olap(3,i)
         id1 = olap(4,i)
         id2 = olap(5,i)

         vx1 = (var0(idg,id1+1) - var0(idg,id1-1)) / 
     1         (x(idg,id1+1) - x(idg,id1-1))
         vx2 = (var0(idg,id2+1) - var0(idg,id2-1)) / 
     1         (x(idg,id2+1) - x(idg,id2-1))

         dx1 = x(idg,id1) - x(irg,irp)
         dx2 = x(idg,id2) - x(irg,irp)
         w1  = 1.0/abs(dx1)**wpower
         w2  = 1.0/abs(dx2)**wpower
         dnm = w1*dx1**2 + w2*dx2**2
         vxp = (var0(idg,id1)-var0(irg,irp))*dx1*w1 +
     1         (var0(idg,id2)-var0(irg,irp))*dx2*w2
         vxp = vxp/dnm

c        Left flux F_1p
         dr  = var0(irg,irp) - var0(idg,id1)
         dl  = 2.0*(-vx1*dx1) - dr
         s   = limiter(dl, dr)
         ul  = var0(idg,id1) + 0.25*s*( (1.0-kkk)*dl + (1.0+kkk)*dr )
         
         dl  = var0(irg,irp) - var0(idg,id1)
         dr  = 2.0*(-vxp*dx1) - dl
         s   = limiter(dl, dr)
         ur  = var0(irg,irp) - 0.25*s*( (1.0-kkk)*dl + (1.0+kkk)*dr )
         
         f1  = fluxfunc(ul, ur)

c        Right flux F_p2
         dr  = var0(idg,id2) - var0(irg,irp)
         dl  = 2.0*(vxp*dx2) - dr
         s   = limiter(dl, dr)
         ul  = var0(irg,irp) + 0.25*s*( (1.0-kkk)*dl + (1.0+kkk)*dr )
         
         dl  = var0(idg,id2) - var0(irg,irp)
         dr  = 2.0*(vx2*dx2) - dl
         s   = limiter(dl, dr)
         ur  = var0(idg,id2) - 0.25*s*( (1.0-kkk)*dl + (1.0+kkk)*dr )
         
         f2  = fluxfunc(ul, ur)

c        Flux at the node p
         fp  = fluxfunc(var0(irg,irp), var0(irg,irp))

         res = (f1 - fp)*dx1*w1 + (f2 - fp)*dx2*w2
         res = 2.0*res/dnm

         if(rkmeth.eq.rk4)then
            var(irg,irp) = var_old(irg,irp) - dt*res
         else if(rkmeth.eq.rktvd)then
            var(irg,irp) = ark(irk)*var_old(irg,irp) + 
     1                     brk(irk)*(var0(irg,irp) - dt*res)
         else
            print*,'Unknown RK scheme'
            stop
         endif
      enddo

      return
      end
c------------------------------------------------------------------------------
c     Update points in overlap region using 3 order meshless method
c------------------------------------------------------------------------------
      subroutine UpdateOverlap_meshless3(irk, dt)
      implicit none
      include 'param.h'
      include 'vars.h'
      include 'common.h'
      integer irk
      real    dt

      integer i, irg, irp, idg, id1, id2
      real    vx1, vx2, vxp, dx1, dx2, dnm, f1, f2, fp, res, dl, dr
      real    fx
      real    u1, u2, s, limiter, w1, w2

      do i=1,nopts
c        Receiving grid and point
         irg = olap(1,i)
         irp = olap(2,i)

c        Donor grid and point
         idg = olap(3,i)
         id1 = olap(4,i)
         id2 = olap(5,i)

         vx1 = (var0(idg,id1+1) - var0(idg,id1-1)) / 
     1         (x(idg,id1+1) - x(idg,id1-1))
         vx2 = (var0(idg,id2+1) - var0(idg,id2-1)) / 
     1         (x(idg,id2+1) - x(idg,id2-1))

         dx1 = x(idg,id1) - x(irg,irp)
         dx2 = x(idg,id2) - x(irg,irp)
         w1  = 1.0/abs(dx1)**wpower
         w2  = 1.0/abs(dx2)**wpower
         dnm = w1*dx1**2 + w1*dx2**2
         vxp = dx2**2 * (var0(idg,id1)-var0(irg,irp)) -
     1         dx1**2 * (var0(idg,id2)-var0(irg,irp))
         vxp = vxp/(dx1*dx2*(dx2-dx1))

         dr  = var0(irg,irp) - var0(idg,id1)
         dl  = 2.0*(-vx1*dx1) - dr
         s   = limiter(dl, dr)
         u1  = var0(idg,id1) + 0.25*s*( (1.0-kkk)*dl + (1.0+kkk)*dr )
         f1  = u1

         dr  = var0(idg,id2) - var0(irg,irp)
         dl  = 2.0*(vxp*dx2) - dr
         s   = limiter(dl, dr)
         u2  = var0(irg,irp) + 0.25*s*( (1.0-kkk)*dl + (1.0+kkk)*dr )
         f2  = u2

         fp  = var0(irg,irp)

         fx  = dx2**2 * (f1 - fp) - dx1**2 * (f2 - fp)
         fx  = fx/(dx1*dx2*(dx2-dx1))
         res = 2.0*fx

         if(rkmeth.eq.rk4)then
            var(irg,irp) = var_old(irg,irp) - dt*res
         elseif(rkmeth.eq.rktvd)then
            var(irg,irp) = ark(irk)*var_old(irg,irp) + 
     1                     brk(irk)*(var0(irg,irp) - dt*res)
         else
            print*,'Unknown RK scheme'
            stop
         endif
      enddo

      return
      end
c------------------------------------------------------------------------------
c     Save solution
c------------------------------------------------------------------------------
      subroutine SaveSol
      implicit none
      include 'param.h'
      include 'vars.h'

      integer i, j

      open(20, file='sol.dat')
      do i=1,ngrid
         do j=1,np(i)
            write(20,*) x(i,j), var(i,j)
         enddo
         write(20,*)
      enddo
      close(20)

      return
      end
c------------------------------------------------------------------------------
c     Plot solution using Grafic
c------------------------------------------------------------------------------
      subroutine PlotGrafic(time)
      implicit none
      include 'param.h'
      include 'vars.h'
      real time

      integer npts, i, j, nper(ngridmax), ilin(ngridmax),
     1        isym(ngridmax), indgr
      real    xx(5*ngridmax*npmax), yy(5*ngridmax*npmax)
      real    InitCond

      npts = 0

      do i=1,ngrid
         do j=1,np(i)
            npts     = npts + 1
            xx(npts) = x(i,j)
            yy(npts) = var(i,j)
         enddo
         nper(i) = np(i)
      enddo

c     Plot error
      do i=1,ngrid
         do j=1,np(i)
            npts     = npts + 1
            xx(npts) = x(i,j)
            yy(npts) = InitCond(x(i,j)-time) - var(i,j)
         enddo
         nper(i+ngrid) = np(i)
      enddo

c     Plot exact solution
      do i=1,1000
         npts = npts + 1
         xx(npts) = (i-1)*0.002
         yy(npts) = InitCond( xx(npts) - time )
      enddo
      nper(5) = 1000

      ilin(1) =  0
      isym(1) = +2

      ilin(2) = 0
      isym(2) = +5

      ilin(3) =  0
      isym(3) = +6

      ilin(4) = 0
      isym(4) = +8

      ilin(5) = 1
      isym(5) = 0

      indgr = 1 + 4 + 16 + 64

#ifdef _GRAFIC
c     call SLEEP(3)
c     call GRSCPT(1,'x')
      call GRLINE(ILIN, ISYM, 5, '~x~u~Linear advection', INDGR, xx, yy,
     1            nper)
#endif

      return
      end
c------------------------------------------------------------------------------
c     Compute maximum and L2 errors
c------------------------------------------------------------------------------
      subroutine Errors(time, errm, err2)
      implicit none
      include 'param.h'
      include 'vars.h'
      real time, errm, err2

      integer i, j, nptot
      real    err
      real    InitCond

      errm = 0.0
      err2 = 0.0
      nptot= 0
      do i=1,ngrid
         do j=1,np(i)
            err = InitCond(x(i,j)-time) - var(i,j)
            err = abs(err)
            errm = max(errm, err)
            err2 = err2 + dx(i)*err**2
            nptot= nptot + 1
         enddo
      enddo

      err2 = sqrt(err2)

      return
      end
