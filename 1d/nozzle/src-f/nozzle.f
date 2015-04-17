      program main
      implicit none
      include 'param.inc'
      integer nf, nc, iter, i, j, nl, nr, pf
      parameter(nf=101)
      parameter(nc=nf-1)
      double precision x(nf), a(nf), q(3,nc), qold(3,nc), res(3,nc),
     +                 dt, residue, residue1, maxres, ar, cost, dum,
     +                 ptarg(nc)

      call init
      print*,'Flux for flow solver =',flux1
      call read_shape(nf, x, a)
c     call init_shape(nf, x, a)

C Initialize flow
      do i=1,nc
         q(1,i) = rinf
         q(2,i) = rinf*uinf
         q(3,i) = pinf/gam1 + 0.5d0*rinf*uinf**2
      enddo

      maxres  = 1.0d-10
      residue = 1.0
      iter    = 0
      open(unit=10, file='res.dat')

C Time step loop
      do while(iter .lt. maxiter .and. residue .gt. maxres)
         call timestep(nc, q, dt)
         call save_old(nc, q, qold)

         do i=1,nc
            do j=1,3
               res(j,i) = 0.0d0
            enddo
         enddo

         call residu(nf, nc, a, q, res)

         do i=1,nc
            ar = 0.5d0*(a(i) + a(i+1))
            do j=1,3
               q(j,i) = qold(j,i) - dt*res(j,i)/dx/ar
            enddo
         enddo

         call sol_residue(nc, res, residue)
         iter = iter + 1
         if(iter .eq. 1) residue1 = residue
         residue = residue/residue1
         write(10,*) iter,dlog10(residue),dt
      enddo
      close(10)

      print*,'Number of flow iterations =',iter

C Read target pressure from a file
      pf = 12
      open(unit=pf, file='flow-target.dat')
      do i=1,nc
         read(pf,*) dum, dum, dum, ptarg(i), dum
      enddo
      close(pf)

C Compute cost function - L2 norm of pressure difference
      cost = 0.0d0
      do i=1,nc
         call costfunc(q(1,i), ptarg(i), cost)
      enddo
      print*,'Cost function =', cost
      open(15, file='cost.dat')
      write(15,*) cost
      close(15)
c     print*, dc, dd, cost

C Save solution into a file
      call result(nf, nc, x, q)

      stop
      end
