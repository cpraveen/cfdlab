      subroutine Solve
      implicit none
      include 'vars.h'

      integer          i, nrks, iter
      double precision time, varmin, varmax

      print*,'Starting solution ...'

      nrks = 3

      time = 0.0d0
      iter = 0

      open(20, file='hist.dat')
      do while(time < Ttime)
         time = time + dt
         iter = iter + 1
         var_old(1:npts) = var(1:npts)
         do i=1,nrks
            call Gradients
            call Residual
            call Update(i)
         enddo

         varmin = 1.0d20
         varmax =-1.0d20
         do i=1,np0
            varmin = dmin1(varmin, var(i))
            varmax = dmax1(varmax, var(i))
         enddo

         write(*,10) iter,time,varmin,varmax
         write(20,10) iter,time,varmin,varmax
10       format(i8,3e16.6)

         if(mod(iter,10).eq.0 .and. makeanim.eq.1)then
            call WriteGNU
            call system('gnuplot cont.gnu')
            call system('gifsicle cont.gif cont0.gif > out.gif')
            call system('mv out.gif cont.gif')
         endif

      enddo
      close(20)

      print*,'Solution range min, max:'
      print*,' Initial=',varmin0,varmax0
      print*,' Final  =',varmin,varmax

      return
      end
