c------------------------------------------------------------------------------
c     1-d spectral element code for heat equation
c     To compile:
c        make
c     To run:
c        ./sem1d
c     Plot results:
c        gnuplot plt.gnu
c        xpdf sol.pdf
c
c     Initial condition: u(0,x) = sin(2*pi*x)
c     Exact solution:    u(t,x) = sin(2*pi*x) * exp(-4*pi*pi*t)
c     Boundary conditions: u(t,0) = u(t,1) = 0
c------------------------------------------------------------------------------
      program main
      include 'SIZE'
      include 'MSH'
      include 'MAT'
      include 'SOL'
      include 'QUAD'
      include 'WORK'

      nx  = lx
      nel = lel

c     create uniform mesh
      xmin = 0.0
      xmax = 1.0
      d    = (xmax - xmin)/nel

      do ie=1,nel
         coord(1,ie) = xmin + (ie-1)*d
         coord(2,ie) = coord(1,ie) + d
         dx(ie) = d
      enddo

c     compute LGL weights and nodes
      n = nx - 1
      call zelegl(n,xg,fg)
      call welegl(n,xg,fg,wt)
      do i=1,nx
         print*,i,xg(i),fg(i),wt(i)
      enddo

c     Compute differentiation matrix
      print*,'Computing differentiation matrix'
      call dmlegl(n,lx-1,xg,fg,dm)

c     compute mass matrix
      call compute_mass_matrix(nx, nel, bm)
      call gsum(nx,nel,bm)

c     set initial condition
      call set_ic(nx, nel, v)

      call output(nx,nel,v)
      call system("cp sol.dat sol0.dat")

      nt = nx*nel

      Tf = 0.01
      t = 0.0
      dt = 1.0e-4
      dt = 0.1*(dx(1)/nx)**2
      do while(t.lt.Tf)
c        dv/dt = A*v
c        res = A*v -------
         call dssum(nx,nel,v)
         call grad(nx,nel,v,vx)
         call compute_element_residual(nx,nel,vx,res)
         call gsum(nx,nel,res)
         res(1) = 0.0
         res(nt) = 0.0
         res(1:nt) = res(1:nt)/bm(1:nt) 
c        res = A*v -------
         v(1:nt) = v(1:nt) + dt*res(1:nt)
         t = t + dt
         print*,t
      enddo

c     Exact solution
      call set_ic(nx, nel, vex)
      call output(nx,nel,vex)
      call system("cp sol.dat ex.dat")

c     numerical solution
      call output(nx,nel,v)

      stop
      end
c------------------------------------------------------------------------------
      subroutine compute_mass_matrix(nx, nel, bm)
      include 'SIZE'
      include 'MSH'
      include 'QUAD'
      real bm(lx,lel)

      print*,'Computing mass matrix ...'

      do ie=1,nel
         do i=1,nx
            bm(i,ie) = 0.5 * dx(ie) * wt(i)
         enddo
      enddo

      return
      end
c------------------------------------------------------------------------------
      subroutine useric()
      include 'MATH'
      include 'WORK'

      u = sin(2.0*pi*x) * exp(-4.0*pi*pi*t)

      return
      end
c------------------------------------------------------------------------------
c     set initial condition by interpolation
      subroutine set_ic(nx, nel, v)
      include 'SIZE'
      include 'MSH'
      include 'QUAD'
      include 'WORK'
      real v(lx,lel)

      do ie=1,nel
         do i=1,nx
            x = coord(1,ie) + 0.5 * dx(ie) * (xg(i) + 1.0)
            call useric()
            v(i,ie) = u
         enddo
      enddo

      return
      end
c------------------------------------------------------------------------------
      subroutine gsum(nx, nel, v)
      include 'SIZE'
      real v(lx,lel)

      do ie=1,nel-1
         atmp = v(nx,ie) + v(1,ie+1)
         v(nx,ie) = atmp
         v(1,ie+1)= atmp
      enddo

      return
      end
c------------------------------------------------------------------------------
c Make v continuous across element boundary by averaging
      subroutine dssum(nx, nel, v)
      include 'SIZE'
      real v(lx,lel)

      do ie=1,nel-1
         atmp = 0.5*(v(nx,ie) + v(1,ie+1))
         v(nx,ie) = atmp
         v(1,ie+1)= atmp
      enddo

      return
      end
c------------------------------------------------------------------------------
      subroutine output(nx,nel,v)
      include 'SIZE'
      include 'QUAD'
      include 'MSH'
      dimension v(lx,lel)

      open(10,file='sol.dat')
      do ie=1,nel
         do i=1,nx
            x = coord(1,ie) + 0.5 * dx(ie) * (xg(i) + 1.0)
            write(10,*) x, v(i,ie)
         enddo
      enddo
      close(10)

      return
      end
c------------------------------------------------------------------------------
      subroutine grad(nx,nel,v,vx)
      include 'SIZE'
      include 'MSH'
      include 'QUAD'
      dimension v(lx,lel), vx(lx,lel)

      do ie=1,nel
         vx(1:nx,ie) = (2.0/dx(ie)) * matmul(dm(1:nx,1:nx) , v(1:nx,ie))
      enddo

      return
      end
c------------------------------------------------------------------------------
      subroutine compute_element_residual(nx,nel,vx,res)
      include 'SIZE'
      include 'QUAD'
      dimension vx(lx,lel), res(lx,lel)

      do ie=1,nel
         do i=1,nx
            res(i,ie) = 0.0
            do j=1,nx
               res(i,ie) = res(i,ie) - vx(j,ie)*dm(j,i)*wt(j)
            enddo
         enddo
      enddo

      return
      end
