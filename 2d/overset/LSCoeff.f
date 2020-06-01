c======================================================================
      subroutine LSCoeff(n,x0,y0,x1,y1,cx,cy)
c======================================================================
c Find coefficients for first derivative using least squares
c n           = Number of neighbours
c x0,y0       = coordinates of point at which grad is required
c x1(n),y1(n) = coordinates of neighbouring points
c cx(n),cy(n) = coefficients
c======================================================================
      implicit none
      integer n
      real    x0, y0, x1(*), y1(*), cx(*), cy(*)

c     local variables
      integer i, j
      real    sdx2, sdy2, sdxdy, dx, dy, dr, w(n), a(n,2), det, 
     +        ainv(2,2)

      sdx2  = 0.0
      sdy2  = 0.0
      sdxdy = 0.0

      do i=1,n
         dx     = x1(i) - x0
         dy     = y1(i) - y0
         dr     = sqrt(dx**2 + dy**2)
         w(i)   = 1.0/dr
         a(i,1) = w(i)*dx
         a(i,2) = w(i)*dy
         sdx2   = sdx2  + a(i,1)**2
         sdy2   = sdy2  + a(i,2)**2
         sdxdy  = sdxdy + a(i,1)*a(i,2)
      enddo

c     invert 2x2 symmetric least squares matrix
      det       = sdx2*sdy2 - sdxdy**2
      ainv(1,1) = sdy2/det
      ainv(1,2) =-sdxdy/det
      ainv(2,1) = ainv(1,2)
      ainv(2,2) = sdx2/det

c     finally the coefficients,
c     [ cx ]
c     [    ] = ainv * a^T * diag(w)
c     [ cy ]
      do i=1,n
         cx(i) = 0.0
         cy(i) = 0.0
         do j=1,2
            cx(i) = cx(i) + ainv(1,j)*a(i,j)*w(i)
            cy(i) = cy(i) + ainv(2,j)*a(i,j)*w(i)
         enddo
      enddo

      return
      end
