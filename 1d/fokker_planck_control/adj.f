c--------------------------------------------------------------------
      program main
      implicit none
      include 'common.inc'
      integer fid, i
      real    dx, x(nc), xx
      real    f1(nc), f2(nc), f1t(nc), f2t(nc), cost
      real    u(2), tstart, tend, tf
      real    ub(2), f1b(nc), f2b(nc), costb

      call setup(tstart, tend, u, dx, x)

      tf   = tend - tstart

      if(initc.eq.0)then
         print*,'Setting initial condition'
         call set_initial_condition(x, f1, f2)
      else
         print*,'Reading initial condition from file'
         fid = 10
         open(fid,file='init.dat',status='old')
         do i=1,nc
            read(fid,*) xx, f1(i), f2(i)
         enddo
         close(fid)
      endif

      call f_target (tend, x, f1t, f2t)

      cost  = 0.0
      costb = 1.0
      f1b   = 0.0
      f2b   = 0.0
      ub    = 0.0
      call COST_FUNCTION_B(u, ub, tf, dx, x, f1, f1b, f2, f2b
     +                           , f1t, f2t, cost, costb)
      print*,"Derivatives = ", ub(1), ub(2)

      fid = 11
      open(fid,file='gradient.dat')
      write(fid,*) ub(1)
      write(fid,*) ub(2)
      close(fid)

      stop
      end
