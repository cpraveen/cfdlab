!------------------------------------------------------------------------------
! Central scheme of Kurganov-Tadmor
! Solve sod problem for 1d Euler equation
!------------------------------------------------------------------------------
module Common
   integer,parameter :: nvar = 3
   integer,parameter :: ng = 3
   real :: gam
   real :: dx
   real :: dt
end module Common

!------------------------------------------------------------------------------
program main
   use Common
   implicit none
   gam = 1.4
   call solve(100, 0.2, 0.5)
end program main

!------------------------------------------------------------------------------
subroutine solve(nc, Tf, cfl)
   use Common
   implicit none
   integer :: nc
   real    :: Tf, cfl
   ! Local variables
   real,parameter :: xmin = 0.0, xmax = 1.0
   integer :: nf, i, iter
   real    :: t, rho, v, p
   real,allocatable :: xc(:), xf(:)
   real,allocatable :: uc(:,:), uf(:,:)
   real,external :: timestep

   nf = nc + 1
   allocate(xc(nc), xf(nf))
   allocate(uc(nvar, -ng+1:nc+ng))
   allocate(uf(nvar, -ng+1:nf+ng))

   ! Make grid
   dx = (xmax - xmin) / nc
   do i=1,nc
      xc(i) = xmin + (i-1)*dx + 0.5*dx
   enddo

   do i=1,nf
      xf(i) = xmin + (i-1)*dx
   enddo

   ! Set initial condition
   do i=1,nf
      if(xf(i) < 0.5)then
         rho = 1.0
         v   = 0.0
         p   = 1.0
      else
         rho = 0.125
         v   = 0.0
         p   = 0.1
      endif
      call prim2con(rho,v,p,uf(:,i))
   enddo

   iter = 0; t = 0.0
   do while(t < Tf)
      ! Update uc
      dt = cfl * timestep(nf, uf)
      dt = min(dt, Tf-t)
      call update(nf, uf, nc, uc)
      t = t + dt
      ! Update uf
      dt = cfl * timestep(nc, uc)
      dt = min(dt, Tf-t)
      call update(nc, uc, nf, uf)
      iter = iter + 1
      t = t + dt
      write(*,*) iter, t
   enddo

   ! Save uf
   open(10, file='sol.txt')
   do i=1,nf
      call con2prim(uf(:,i), rho, v, p)
      write(10,*)xf(i), rho, v, p
   enddo
   close(10)
end subroutine solve

!------------------------------------------------------------------------------
real function timestep(n, u)
   use Common
   implicit none
   integer :: n
   real :: u(nvar, -ng+1:n+ng)
   ! Local variables
   integer :: i
   real :: rho, v, p

   timestep = 1.0e20
   do i=1,n
      call con2prim(u(:,i), rho, v, p)
      timestep = min(timestep, dx/(abs(v) + sqrt(gam*p/rho)))
   enddo

end function timestep

!------------------------------------------------------------------------------
subroutine prim2con(rho, v, p, u)
   use Common
   implicit none
   real,intent(in) :: rho, v, p
   real,intent(inout) :: u(nvar)
   u(1) = rho
   u(2) = rho * v
   u(3) = p/(gam - 1.0) + 0.5 * rho * v**2
end subroutine prim2con

!------------------------------------------------------------------------------
subroutine con2prim(u, rho, v, p)
   use Common
   implicit none
   real,intent(in) :: u(nvar)
   real,intent(inout) :: rho, v, p
   rho = u(1)
   v   = u(2) / u(1)
   p   = (gam - 1.0) * (u(3) - 0.5 * u(2)**2 / u(1))
end subroutine con2prim

!------------------------------------------------------------------------------
subroutine jacobian(rho, v, p, A)
   use Common
   implicit none
   real,intent(in) :: rho, v, p
   real,intent(inout) :: A(nvar,nvar)
   ! Local variables
   real :: v2, h

   v2 = v**2
   h = gam * p / (rho * (gam - 1)) + 0.5 * v2

   A(1,1) = 0.0
   A(2,1) = 0.5 * (gam - 3.0) * v2
   A(3,1) = (0.5 * (gam - 1.0) * v2 - h) * v

   A(1,2) = 1.0
   A(2,2) = (3.0 - gam) * v
   A(3,2) = h - (gam - 1.0) * v2

   A(1,3) = 0.0
   A(2,3) = gam - 1.0
   A(3,3) = gam * v
end subroutine jacobian

!------------------------------------------------------------------------------
real function minmod(ujm1, uj, ujp1)
   implicit none
   real :: ujm1, uj, ujp1
   ! Local variables
   real :: db, df, sb, sf
   db = uj - ujm1
   df = ujp1 - uj
   sb = sign(1.0, db)
   sf = sign(1.0, df)
   if(sb == sf)then
      minmod = sb * min(abs(db), abs(df))
   else
      minmod = 0.0
   endif
end function minmod

!------------------------------------------------------------------------------
! Given u1, update u2
subroutine update(n1, u1, n2, u2)
   use Common
   implicit none
   integer :: n1, n2
   real :: u1(nvar, -ng+1:n1+ng)
   real :: u2(nvar, -ng+1:n2+ng)
   ! Local variables
   integer :: i, j, il, ir, l, r
   real :: rho, v, p, lam
   real :: ux(nvar, 0:n1+1), f(nvar, 0:n1+1), uh(nvar), fx(nvar), A(nvar,nvar)
   real,external :: minmod

   ! Fill ghost values: neumann bc
   do i=1,ng
      u1(:,1-i ) = u1(:,1 )
      u1(:,n1+i) = u1(:,n1)
   enddo

   do i=1,n1
      do j=1,nvar
         ux(j,i) = minmod(u1(j,i-1), u1(j,i), u1(j,i+1)) / dx
      enddo
      call con2prim(u1(:,i), rho, v, p)
      call jacobian(rho, v, p, A)
      fx = matmul(A, ux(:,i))
      uh = u1(:,i) - 0.5 * dt * fx
      call con2prim(uh, rho, v, p)
      f(1,i) = rho * v
      f(2,i) = p + rho * v**2
      f(3,i) = (gam * p / (gam - 1.0) + 0.5 * rho * v**2) * v
   enddo

   ! Fill ghost values: neumann bc
   do i=1,ng
      ux(:,0 )   = ux(:,1 )
      ux(:,n1+1) = ux(:,n1)
      f(:,0 )   = f(:,1 )
      f(:,n1+1) = f(:,n1)
   enddo

   lam = dt/dx
   if(n2 > n1)then ! update uf
      il = -1
      ir =  0
   else            ! update uc
      il =  0
      ir = +1
   endif

   do i=1,n2
      l = i + il
      r = i + ir
      do j=1,nvar
         u2(j,i) = 0.5*(u1(j,l) + u1(j,r)) &
                   + (dx/8.0) * (ux(j,l) - ux(j,r)) &
                   - lam * (f(j,r) - f(j,l))
      enddo
   enddo
end subroutine update
