      program main
      implicit none
      include 'param.inc'
      integer nf, nc, iter, i, j, nl, nr, pf, irk
      integer nfmax, ncmax
      parameter(nfmax=1001)
      parameter(ncmax=nfmax-1)
      double precision x(nfmax), a(nfmax), q(3,ncmax), qold(3,ncmax), 
     +                 res(3,ncmax),
     +                 dt, residue, residue1, maxres, ar, cost, dum,
     +                 ptarg(ncmax), p, xx

      call init(nf, nc)
      print*,'Flux for flow solver =',flux1
c     call read_shape(nf, x, a)
      call init_shape(nf, x, a)

C Initialize flow
      do i=1,nc
         q(1,i) = rinf
         q(2,i) = rinf*uinf
         xx = 0.5d0*(x(i)+x(i+1))
         p  = pinf + (xx-xmin)*(pout-pinf)/L
         q(3,i) = p/gam1 + 0.5d0*rinf*uinf**2
      enddo

      maxres  = 1.0d-10
      residue = 1.0
      iter    = 0
      open(unit=10, file='res.dat')

C Time step loop
      do while(iter .lt. maxiter .and. residue .gt. maxres)
         call timestep(nc, q, dt)
         call save_old(nc, q, qold)

         do irk=1,3

         do i=1,nc
            do j=1,3
               res(j,i) = 0.0d0
            enddo
         enddo

         call residu(nf, nc, a, q, res)

         do i=1,nc
            ar = 0.5d0*(a(i) + a(i+1))
            do j=1,3
               q(j,i) = ark(irk)*qold(j,i) +
     1                  (1.0d0-ark(irk))*(q(j,i) - dt*res(j,i)/dx/ar)
            enddo
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
c     pf = 12
c     open(unit=pf, file='flow-target.dat')
c     do i=1,nc
c        read(pf,*) dum, dum, dum, ptarg(i), dum
c     enddo
c     close(pf)

C Compute cost function - L2 norm of pressure difference
      call costfunc(nc, q, cost)
      print*,'Cost function =', cost
      open(15, file='cost.dat')
      write(15,'(e24.15)') cost
      close(15)
c     print*, dc, dd, cost

C Save solution into a file
      call result(nf, nc, x, q)

      stop
      end
