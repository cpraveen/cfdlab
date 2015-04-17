C Computes finite volume residual
      subroutine residu(nf, nc, a, q, res)
      implicit none
      include 'param.inc'
      integer nf, nc
      double precision a(nf), q(3,nc), res(3,nc)

      integer i, j, nl, nr, pf

      call flux_in(a(1), res(1,1))
      do i=2,nf-1
         nl = i-1
         nr = i
         if(flux1 .eq. 1)then
            call flux_ausm(a(i), q(1,nl), q(1,nr), res(1,nl), res(1,nr))
         elseif(flux1 .eq. 2)then
            call flux_kfvs(a(i), q(1,nl), q(1,nr), res(1,nl), res(1,nr))
         elseif(flux1 .eq. 3)then
            call flux_lf(a(i), q(1,nl), q(1,nr), res(1,nl), res(1,nr))
         endif
      enddo
      call flux_out(a(nf), q(1,nf-1), q(1,nf-2), res(1,nf-1))

      do i=1,nc
         call source_term(a(i), a(i+1), q(1,i), res(1,i))
      enddo

      return
      end
