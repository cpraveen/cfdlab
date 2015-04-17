      program main
      implicit none
c     include 'param.inc'
      include 'param_bmach.inc'
      integer nf, nc, iter, i, j, nl, nll, nr, nrr, pf, irk
      integer nfmax, ncmax
      parameter(nfmax=1001)
      parameter(ncmax=nfmax-1)
      double precision x(nfmax), a(nfmax), q(3,ncmax), qold(3,ncmax), 
     +                 res(3,ncmax),
     +                 dt, residue, residue1, maxres, ar, cost, dum,
     +                 ptarg(ncmax)
      double precision qb(3,ncmax), resb(3,ncmax), resbold(3,ncmax),
     +                 ab(nfmax), 
     +                 qb1(3,ncmax), costb
      double precision dJdc, dJdd, dcb, ddb, nozzleb
      double precision ql(3), qr(3), qlb(3), qrb(3), qinfb(3)
      double precision resbl(3), resbr(3)

      call init(nf, nc)
      print*,'Flux for adjoint solver =',flux2
      call read_flow(nc, q)
c     call read_shape(nf, x, a)
      call init_shape(nf, x, a)
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
c     print*,'Reading target pressure'
c     pf = 12
c     open(unit=pf, file='flow-target.dat')
c     do i=1,nc
c        read(pf,*) dum, dum, dum, ptarg(i), dum
c     enddo
c     close(pf)

C Compute cost function - L2 norm of pressure difference
      costb = 1.0d0
      call costfunc_bq(nc, q, qb1, cost, costb)

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

         do irk=1,3

         do i=1,nc
            qb(1,i) = 0.0d0
            qb(2,i) = 0.0d0
            qb(3,i) = 0.0d0
         enddo

         do i=1,nc
            call source_term_bq(a(i), a(i+1), q(1,i), qb(1,i), 
     +                          res(1,i), resb(1,i))
         enddo

c        call flux_out_bq(a(nf), q(1,nf-1), qb(1,nf-1),
c    +                    q(1,nf-2), qb(1,nf-2), res(1,nf-1),
c    +                    resb(1,nf-1))
         call flux_out2_bq(a(nf), q(1,nf-1), qb(1,nf-1),
     +                    res(1,nf-1),
     +                    resb(1,nf-1))
c     do i=1,nf-1
      do i=2,nf-1
         nl = i-1
         nr = i
         nll= i-2
         nrr= i+1
         if(i.eq.1)then
            call muscl(qinf    , qinf   , q(1,nr), ql)
            call muscl(q(1,nrr), q(1,nr), qinf, qr)
         elseif(i.eq.2)then
            ql(:) = q(:,nl)
c           call muscl(qinf    , q(1,nl), q(1,nr), ql)
            call muscl(q(1,nrr), q(1,nr), q(1,nl), qr)
         elseif(i.eq.nf-1)then
            ql(:) = q(:,nl)
            qr(:) = q(:,nr)
         else
            call muscl(q(1,nll), q(1,nl), q(1,nr), ql)
            call muscl(q(1,nrr), q(1,nr), q(1,nl), qr)
         endif

            qlb(:) = 0.0d0
            qrb(:) = 0.0d0
            if(flux2 .eq. 1)then
            call flux_ausm_bq(a(i), ql, qlb, qr, qrb,
     +                   res(1,nl), resb(1,nl), res(1,nr), resb(1,nr))
            elseif(flux2 .eq. 2)then
            call flux_kfvs_bq(a(i), ql, qlb, qr, qrb,
     +                   res(1,nl), resb(1,nl), res(1,nr), resb(1,nr))
            elseif(flux2 .eq. 3)then
            call flux_lf_bq(a(i), ql, qlb, qr, qrb,
     +                   res(1,nl), resb(1,nl), res(1,nr), resb(1,nr))
            elseif(flux2 .eq. 4)then
            call flux_roe_bq(a(i), ql, qlb, qr, qrb,
     +                   res(1,nl), resb(1,nl), res(1,nr), resb(1,nr))
            endif

         if(i.eq.1)then
            call muscl_bq(qinf    , qinfb    , qinf   , qinfb, 
     1                    q(1,nr), qb(1,nr), ql, qlb)
            call muscl_bq(q(1,nrr), qb(1,nrr), q(1,nr), qb(1,nr),
     1                    qinf   , qinfb   , qr, qrb)
         elseif(i.eq.2)then
            qb(:,nl) = qb(:,nl) + qlb(:)
c           call muscl_bq(qinf    , qinfb    , q(1,nl), qb(1,nl), 
c    1                    q(1,nr), qb(1,nr), ql, qlb)
            call muscl_bq(q(1,nrr), qb(1,nrr), q(1,nr), qb(1,nr),
     1                    q(1,nl), qb(1,nl), qr, qrb)
         elseif(i.eq.nf-1)then
            qb(:,nl) = qb(:,nl) + qlb(:)
            qb(:,nr) = qb(:,nr) + qrb(:)
         else
            call muscl_bq(q(1,nll), qb(1,nll), q(1,nl), qb(1,nl), 
     1                    q(1,nr), qb(1,nr), ql, qlb)
            call muscl_bq(q(1,nrr), qb(1,nrr), q(1,nr), qb(1,nr),
     1                    q(1,nl), qb(1,nl), qr, qrb)
         endif
      enddo
      call flux_in_steger_bq(a(1), q(1,1), qb(1,1), res(1,1), resb(1,1))

         do i=1,nc
            ar = 0.5d0*(a(i) + a(i+1))
            do j=1,3
               qb(j,i)   = qb(j,i) + qb1(j,i)
c              resb(j,i) = resbold(j,i) - dt*qb(j,i)/dx/ar
               resb(j,i) = ark(irk)*resbold(j,i) +
     1                   (1.0d0-ark(irk))*(resb(j,i) - dt*qb(j,i)/dx/ar)
            enddo
         enddo

         enddo

         call sol_residue(nc, qb, residue)
         iter = iter + 1
         if(iter .eq. 1) residue1 = residue
         write(10,*) iter, dlog10(residue)
      enddo
      close(10)

      print*,'Number of adjoint iterations =',iter

c derivative wrt inflow mach number
      call flux_in_steger_bmach(a(1), q(1,1),
     1                          res(1,1), resb(1,1))
      print*,'Derive wrt minf =', minfb


C Save adjoint solution into a file
      ab(1)  = 0.0d0
      ab(nf) = 0.0d0
      call result_adj(nf, nc, x, ab, resb)
      print*,'Adjoint solution written into flowb.dat and shapeb.dat'

      stop
      end
