C Computes finite volume residual
      subroutine residu(nf, nc, a, q, res)
      implicit none
      include 'param.inc'
      integer nf, nc
      double precision a(*), q(3,*), res(3,*)

      integer i, j, nl, nll, nr, nrr, pf
      double precision ql(3), qr(3)

c     call flux_in(a(1), res(1,1))
c     call flux_in2(a(1), q(1,1), res(1,1))
      call flux_in_steger(a(1), q(1,1), res(1,1))
      do i=2,nf-1
c     do i=1,nf-1
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
         if(flux1 .eq. 1)then
            call flux_ausm(a(i), ql, qr, res(1,nl), res(1,nr))
         elseif(flux1 .eq. 2)then
            call flux_kfvs(a(i), ql, qr, res(1,nl), res(1,nr))
         elseif(flux1 .eq. 3)then
            call flux_lf(a(i), ql, qr, res(1,nl), res(1,nr))
         elseif(flux1 .eq. 4)then
            call flux_roe(a(i), ql, qr, res(1,nl), res(1,nr))
         endif
      enddo
      call flux_out2(a(nf), q(1,nf-1), res(1,nf-1))

      do i=1,nc
         call source_term(a(i), a(i+1), q(1,i), res(1,i))
      enddo

      return
      end
