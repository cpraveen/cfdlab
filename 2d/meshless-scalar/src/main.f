c     2-d meshless solver for linear advection equation with periodic
c     boundary conditions
      program main
      implicit none
      include 'vars.h'

c     ===========set some parameters===========
c     domain
      xmin = 0.0d0
      xmax = 1.0d0
      ymin = 0.0d0
      ymax = 1.0d0

c     grid size
      nx   = 101
      ny   = 101

c     perturb the grid, 0=no, 1=yes
      perturb = 0
      pertx   = 0.1d0
      perty   = 0.1d0

c     1=first order, 2=higher, 3=higher+limiter
      order= 2

c     factor in muscl reconstruction
      kkk = 1.0d0/3.0d0

c     initial condition, see SetInitCond.f
      ictype = 1

c     advection velocity, see AdvecVel.f
      veltype = 1

c     time step
      dt = 0.005d0

c     final time
      Ttime = 1.0d0

c     animation req gnuplot and gifsicle, set makeanim=0 to diable this
      makeanim = 1
      animinterval = 10

c     ===========solution begins from here===========
c     generate grid and connectivity
      call GenPtsConn

c     compute coefficients for deriv using least-squares
      call LeastSquares

c     initialize solution
      call SetInitCond
      if(makeanim.eq.1)then
         call WriteGNU
         call system('gnuplot cont.gnu')
         call system('mv cont0.gif cont.gif')
      endif

c     find global time-step
      call InitTimeStep

c     solve
      call Solve

c     save final solution
      call WriteResult

c     remove some files
      call system('rm -f cont0.gif cnt.dat')

      print*,'Output files: hist.dat result.dat result.vig cont.gif'

      stop
      end
