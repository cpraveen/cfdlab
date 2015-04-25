c--------------------------------------------------------------------
c Solve linear convection equation with periodic bc
c Author: Praveen. C, http://praveen.tifrbng.res.in
c
c Test case from
c Annunziato & Borzi: Euro. Jnl. of Applied Mathematics (2014)
c vol. 25, pp. 1-15
c
c     f1_t + (a1(x) f1)_x = -mu*f1 + mu*f2
c     f2_t + (a2(x) f2)_x =  mu*f1 - mu*f2
c
c  mu is set in common.inc
c  a1(x) and a2(x) can be set in the functions below
c  For initial condition, edit subroutine set_initial_condition
c  Solution is saved into file sol.dat; plot in gnuplot using
c     bash$ gnuplot plot.gnu
c     bash$ python plot.py
c     bash$ xpdf sol.pdf
c--------------------------------------------------------------------
      program main
      implicit none
      include 'common.inc'
      integer fid, i
      real    dx, x(nc), xx
      real    f1(nc), f2(nc), f1t(nc), f2t(nc), cost
      real    u(2), tstart, tend, tf

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
      call output(x, f1t, f2t)
      call system("cp sol.dat target.dat")

      cost = 0.0
      call cost_function(u, tf, dx, x, f1, f2, f1t, f2t, cost)

      fid = 20
      open(fid,file='cost.dat')
      write(fid,*) cost
      close(fid)

      stop
      end
