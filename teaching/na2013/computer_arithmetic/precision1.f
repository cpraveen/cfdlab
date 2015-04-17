c Even a simple number like 0.1 cannot be exactly represented in the 
c computer
      implicit none
      real xr
      double precision xd

      xr = 0.1
      xd = 0.1d0

      write(*,'(e26.20)') xr
      write(*,'(e26.20)') xd

      stop
      end
