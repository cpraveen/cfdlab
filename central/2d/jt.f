      program main
      parameter(nx=400,ny=400,mn=4,md=3)
      real u(nx+2*md,ny+2*md,mn)

      xmin = 0.0
      xmax = 1.0
      ymin = xmin
      ymax = xmax
      cfl = 0.475
      gamma = 1.4
      theta = 2.0
      tf = 0.3
      dx = (xmax - xmin)/nx
      dy = (ymax - ymin)/ny

c     Set initial condition
      do j=md+1,md+ny
         do i=md+1,md+nx
            x = xmin + (i-1-md)*dx + 0.5*dx
            y = ymin + (j-1-md)*dy + 0.5*dy
c           Upper, right
            if(x .ge. 0.5 .and. y .ge. 0.5)then 
               rho = 1.5
               vex = 0.0
               vey = 0.0
               pre = 1.5
c           Upper, left
            else if(x .le. 0.5 .and. y .ge. 0.5)then
               rho = 0.5323
               vex = 1.206
               vey = 0.0
               pre = 0.3
c           Lower, left
            else if(x .le. 0.5 .and. y .le. 0.5)then
               rho = 0.138
               vex = 1.206
               vey = 1.206
               pre = 0.029
c           Lower, right
            else
               rho = 0.5323
               vex = 0.0
               vey = 1.206
               pre = 0.3
            endif
            u(i,j,1) = rho
            u(i,j,2) = rho * vex
            u(i,j,3) = rho * vey
            u(i,j,4) = pre/(gamma-1.0) + 0.5*rho*(vex**2 + vey**2)
         enddo
      enddo

      call euler2d(nx,ny,dx,dy,cfl,gamma,theta,tf,u)

      open(10,file='sol.plt')
      write(10,*)"TITLE = Central_scheme"
      write(10,*)"VARIABLES = x, y, density, vex, vey, pre"
      write(10,*)"ZONE STRANDID=1,SOLUTIONTIME=",tf,",I=",nx,",J=",ny,
     &           ",DATAPACKING=POINT"
      do j=md+1,md+ny
         do i=md+1,md+nx
            x = xmin + (i-1-md)*dx + 0.5*dx
            y = ymin + (j-1-md)*dy + 0.5*dy
            rho = u(i,j,1)
            vex = u(i,j,2) / rho
            vey = u(i,j,3) / rho
            pre = (gamma-1.0)*(u(i,j,4) - 0.5*rho*(vex**2 + vey**2))
            write(10,*)x,y,rho,vex,vey,pre
         enddo
      enddo
      close(10)

      end program main

      subroutine EULER2D(nx,ny,dx,dy,cfl,gamma,theta,tf,u)
******************************************************************
* INPUT nx, ny: # of cells in x-, y-direction
* dx, dy: step sizes in x-, y-direction
* cfl: CFL # gamma: adiabatic constant of gas
* tf: final time
* theta=1: MM1 limiter; =2: MM2 limiter; >2: UNO limiter.
* u: initial cell averages of conservative variables.
* Supply entries of u((md+1):(nx+md),(md+1):(ny+md),4)
* OUTPUT u: cell averages at final time "tf"
* REMARK 1. Reset "nxd","nyd" to adjust array dimensions.
* 2. Padded to each side of the computational domain are
* "md" ghost cells, average values on which are
* assigned by boundary conditions.
* Taken from:
*     Jiang and Tadmar: Nonoscillatory central schemes for
*     multidimensional hyperbolic conservation laws
******************************************************************
      parameter(md=3, mn=4)
      real u(nx+2*md,ny+2*md,mn)
      real ux(nx+2*md,ny+2*md,mn), uy(nx+2*md,ny+2*md,mn)
      real f(nx+2*md,ny+2*md,mn), fx(nx+2*md,ny+2*md,mn)
      real g(nx+2*md,ny+2*md,mn), gy(nx+2*md,ny+2*md,mn)
      real v(nx+2*md,ny+2*md), du(nx+2*md,2), df(nx+2*md,2)

      xmin(a,b) = 0.5*(sign(1.,a)+sign(1.,b))*min(abs(a),abs(b))
      xmic(z,a,b) = xmin( z*xmin(a,b), 0.5*(a+b) )

      gm1 = gamma - 1.0
      tc = 0.0
      istop = 0

      do 1000 nt = 1, 1000000
      do 999 io = 0, 1

* Periodic boundary condition in both x- & y-direction
      do 101 m = 1, mn
      do 100 i = 1, md
c     do 100 j = md + 1, ny + md
      do 100 j = 1, ny + 2*md
      u(i, j,m) = u(nx+i,j,m)
      u(nx + md+i,j,m) = u(md+i,j,m)
100   continue
      do 101 j = 1, md
      do 101 i = 1, nx + 2*md
      u(i,j, m) = u(i,ny+j,m)
      u(i,ny + md+j,m) = u(i,md+j,m)
101   continue

* Compute f & g and maximum wave speeds "em_x, em_y".
* See (2.1) & (4.3).
      em_x = 1.e-15
      em_y = 1.e-15
      do 200 j = 1, ny + 2*md
      do 200 i = 1, nx + 2*md
      den = u(i,j,1)
      vex = u(i,j,2) / den
      vey = u(i,j,3) / den
      eng = u(i,j,4)
      pre = gm1 * ( eng - .5*den*( vex*vex + vey*vey ) )
      cvel = sqrt( gamma * pre / den )
      em_x = max( em_x, abs(vex) + cvel )
      em_y = max( em_y, abs(vey) + cvel )
      f(i,j,1) = den * vex
      f(i,j,2) = den * vex**2 + pre
      f(i,j,3) = den * vex * vey
      f(i,j,4) = vex * ( pre + eng )
      g(i,j,1) = den * vey
      g(i,j,2) = den * vex * vey
      g(i,j,3) = den * vey**2 + pre
      g(i,j,4) = vey * ( pre + eng )
200   continue

* Compute numerical derivatives "ux", "uy", "fx", "gy".
* See (3.1) & (4.1)
      do 330 m = 1, mn
      do 310 j = 3, ny + 2*md - 2
      do 301 i = 1, nx + 2*md - 1
      du(i,1) = u(i+1,j,m) - u(i,j,m)
      df(i,1) = f(i+1,j,m) - f(i,j,m)
301   continue
      do 302 i = 1, nx + 2*md - 2
      du(i,2) = du(i+1,1) - du(i,1)
      df(i,2) = df(i+1,1) - df(i,1)
302   continue
      if( theta .lt. 2.5 ) then
      do 303 i = 3, nx + 2*md - 2
      ux(i,j,m) = xmic( theta, du(i-1,1), du(i,1) )
      fx(i,j,m) = xmic( theta, df(i-1,1), df(i,1) )
303   continue
      else
      do 304 i = 3, nx + 2*md - 2
      ux(i,j,m)=xmin(du(i-1,1)+.5*xmin(du(i-2,2),du(i-1,2)),
     &     du(i, 1)-.5*xmin(du(i-1,2),du(i, 2)))
      fx(i,j,m)=xmin(df(i-1,1)+.5*xmin(df(i-2,2),df(i-1,2)),
     &     df(i, 1)-.5*xmin(df(i-1,2),df(i, 2)))
304   continue
      endif
310   continue
      do 320 i = 3, nx + 2*md - 2
      do 311 j = 1, ny + 2*md - 1
      du(j,1) = u(i,j+1,m) - u(i,j,m)
      df(j,1) = g(i,j+1,m) - g(i,j,m)
311   continue
      do 312 j = 1, ny + 2*md - 2
      du(j,2) = du(j+1,1) - du(j,1)
      df(j,2) = df(j+1,1) - df(j,1)
312   continue
      if( theta .lt. 2.5 ) then
      do 313 j = 3, ny + 2*md - 2
      uy(i,j,m) = xmic( theta, du(j-1,1), du(j,1) )
      gy(i,j,m) = xmic( theta, df(j-1,1), df(j,1) )
313   continue
      else
      do 314 j = 3, ny + 2*md - 2
      uy(i,j,m)=xmin(du(j-1,1)+.5*xmin(du(j-2,2),du(j-1,2)),
     & du(j, 1)-.5*xmin(du(j-1,2),du(j, 2)))
      gy(i,j,m)=xmin(df(j-1,1)+.5*xmin(df(j-2,2),df(j-1,2)),
     & df(j, 1)-.5*xmin(df(j-1,2),df(j, 2)))
314   continue
      endif
320   continue
330   continue

* Compute time step size according to the input CFL #.
      if(io.eq.0) then
      dt = cfl / max( em_x/dx, em_y/dy )
      if( ( tc + 2.*dt ) .ge. tf ) then
      dt = 0.5 * ( tf - tc )
      istop = 1
      endif
      endif
      dtcdx2 = 0.5 * dt / dx
      dtcdy2 = 0.5 * dt / dy
* Compute the flux values of f & g at half time step.
* See (2.15) & (2.16).
      do 400 j = 3, ny + 2*md - 2
      do 400 i = 3, nx + 2*md - 2
      den = u(i,j,1) - dtcdx2*fx(i,j,1) - dtcdy2*gy(i,j,1)
      xmt = u(i,j,2) - dtcdx2*fx(i,j,2) - dtcdy2*gy(i,j,2)
      ymt = u(i,j,3) - dtcdx2*fx(i,j,3) - dtcdy2*gy(i,j,3)
      eng = u(i,j,4) - dtcdx2*fx(i,j,4) - dtcdy2*gy(i,j,4)
      pre = gm1 * ( eng - .5 * ( xmt*xmt + ymt*ymt ) / den )
      f(i,j,1) = xmt
      f(i,j,2) = xmt * xmt / den + pre
      f(i,j,3) = xmt * ymt / den
      f(i,j,4) = xmt / den * ( pre + eng )
      g(i,j,1) = ymt
      g(i,j,2) = xmt * ymt / den
      g(i,j,3) = ymt * ymt / den + pre
      g(i,j,4) = ymt / den * ( pre + eng )
400   continue
* Compute the values of "u" at the next time level. See (2.16).
      do 510 m = 1, mn
      do 501 j = md + 1 - io, ny + md - io
      do 501 i = md + 1 - io, nx + md - io
      v(i,j) = 0.25 * ( u(i,j, m) + u(i+1,j, m)
     & + u(i,j+1,m) + u(i+1,j+1,m) )
     & + 0.0625 * ( ux(i,j, m) - ux(i+1,j, m)
     & + ux(i,j+1,m) - ux(i+1,j+1,m)
     & + uy(i,j, m) + uy(i+1,j, m)
     & - uy(i,j+1,m) - uy(i+1,j+1,m) )
     & + dtcdx2 * ( f(i,j, m) - f(i+1,j, m)
     & + f(i,j+1,m) - f(i+1,j+1,m) )
     & + dtcdy2 * ( g(i,j, m) + g(i+1,j, m)
     & - g(i,j+1,m) - g(i+1,j+1,m) )
501   continue
      do 502 j = md + 1, ny + md
      do 502 i = md + 1, nx + md
      u(i,j,m) = v(i-io,j-io)
502   continue
510   continue
      tc = tc + dt
999   continue
      print*,'time =',tc,'dt =',dt
      if(istop.eq.1 ) goto 1001
1000  continue
1001  return
      end
