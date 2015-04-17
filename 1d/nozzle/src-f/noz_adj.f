      program main
      implicit none
      include 'param.inc'
      integer nf, nc, iter, i, j, nl, nr, pf
      parameter(nf=101)
      parameter(nc=nf-1)
      double precision x(nf), a(nf), q(3,nc), qold(3,nc), res(3,nc),
     +                 dt, residue, residue1, maxres, ar, cost, dum,
     +                 ptarg(nc)
      double precision qb(3,nc), resb(3,nc), resbold(3,nc), ab(nf), 
     +                 qb1(3,nc), costb
      double precision dJdc, dJdd, dcb, ddb, nozzleb

      call init
      print*,'Flux for adjoint solver =',flux2
      call read_flow(nc, q)
      call read_shape(nf, x, a)
c     call init_shape(nf, x, a)
      call timestep(nc, q, dt)
      print*,'Time step =',dt

C Initialize flow
      do i=1,nc
         qb1(1,i) = 0.0d0
         qb1(2,i) = 0.0d0
         qb1(3,i) = 0.0d0
         res(1,i) = 0.0d0
         res(2,i) = 0.0d0
         res(3,i) = 0.0d0
      enddo

C Read target pressure from a file
      print*,'Reading target pressure'
      pf = 12
      open(unit=pf, file='flow-target.dat')
      do i=1,nc
         read(pf,*) dum, dum, dum, ptarg(i), dum
      enddo
      close(pf)

C Compute cost function - L2 norm of pressure difference
      costb = 1.0d0
      do i=1,nc
         call costfunc_bq(q(1,i), qb1(1,i), ptarg(i), cost, costb)
      enddo

      maxres  = 1.0d-10
      residue = 1.0
      iter    = 0

      do i=1,nc
         do j=1,3
            resb(j,i) = 0.0d0
         enddo
      enddo

      print*,'Beginning of adjoint iterations'
      open(unit=10, file='resb.dat')
C Time step loop
      do while(iter .lt. maxiter .and. residue .gt. maxres)
         call save_old(nc, resb, resbold)

         do i=1,nc
            qb(1,i) = 0.0d0
            qb(2,i) = 0.0d0
            qb(3,i) = 0.0d0
         enddo

         do i=1,nc
            call source_term_bq(a(i), a(i+1), q(1,i), qb(1,i), 
     +                          res(1,i), resb(1,i))
         enddo

         call flux_out_bq(a(nf), q(1,nf-1), qb(1,nf-1),
     +                    q(1,nf-2), qb(1,nf-2), res(1,nf-1),
     +                    resb(1,nf-1))
         do i=nf-1,2,-1
            nl = i-1
            nr = i
            if(flux2 .eq. 1)then
            call flux_ausm_bq(a(i), q(1,nl), qb(1,nl), q(1,nr), 
     +                   qb(1,nr), res(1,nl), resb(1,nl), res(1,nr),
     +                   resb(1,nr))
            elseif(flux2 .eq. 2)then
            call flux_kfvs_bq(a(i), q(1,nl), qb(1,nl), q(1,nr), 
     +                   qb(1,nr), res(1,nl), resb(1,nl), res(1,nr),
     +                   resb(1,nr))
            elseif(flux2 .eq. 3)then
            call flux_lf_bq(a(i), q(1,nl), qb(1,nl), q(1,nr), 
     +                   qb(1,nr), res(1,nl), resb(1,nl), res(1,nr),
     +                   resb(1,nr))
            endif
         enddo
         call flux_in_bq(a(1), res(1,1), resb(1,1))

         do i=1,nc
            ar = 0.5d0*(a(i) + a(i+1))
            do j=1,3
               qb(j,i)   = qb(j,i) + qb1(j,i)
               resb(j,i) = resbold(j,i) - dt*qb(j,i)/dx/ar
            enddo
         enddo

         call sol_residue(nc, qb, residue)
         iter = iter + 1
         if(iter .eq. 1) residue1 = residue
         write(10,*) iter, dlog10(residue)
      enddo
      close(10)

      print*,'Number of adjoint iterations =',iter

C Computation of contribution of shape perturbation
         do i=1,nf
            ab(i) = 0.0d0
         enddo

         do i=1,nc
            call source_term_ba(a(i), ab(i), a(i+1), ab(i+1), q(1,i), 
     +                          res(1,i), resb(1,i))
         enddo

         call flux_out_ba(a(nf), ab(nf), q(1,nf-1), q(1,nf-2), 
     +                    res(1,nf-1), resb(1,nf-1))
         do i=nf-1,2,-1
            nl = i-1
            nr = i
            if(flux2 .eq. 1)then
            call flux_ausm_ba(a(i), ab(i), q(1,nl), q(1,nr), 
     +                   res(1,nl), resb(1,nl), res(1,nr), resb(1,nr))
            elseif(flux2 .eq. 2)then
            call flux_kfvs_ba(a(i), ab(i), q(1,nl), q(1,nr), 
     +                   res(1,nl), resb(1,nl), res(1,nr), resb(1,nr))
            elseif(flux2 .eq. 3)then
            call flux_lf_ba(a(i), ab(i), q(1,nl), q(1,nr), 
     +                   res(1,nl), resb(1,nl), res(1,nr), resb(1,nr))
            endif
         enddo
         call flux_in_ba(a(1), ab(1), res(1,1), resb(1,1))

C Save adjoint solution into a file
      ab(1)  = 0.0d0
      ab(nf) = 0.0d0
      call result_adj(nf, nc, x, ab, resb)
      print*,'Adjoint solution written into flowb.dat and shapeb.dat'

      stop

C Run this part if you want gradients with respect to the parameters c
C and d in the definition of nozzle shape
      dJdc = 0.0d0
      dJdd = 0.0d0
      do i=1,nf
         nozzleb = 1.0d0
         call NOZZLE_BA(L, ain, aout, x(i), dc, dcb, dd, ddb, nozzleb)
         dJdc = dJdc + ab(i)*dcb
         dJdd = dJdd + ab(i)*ddb
      enddo
      print*,'dJdc =',dJdc
      print*,'dJdd =',dJdd

      stop
      end
