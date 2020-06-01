c======================================================================
c Generate points and connectivity
      subroutine GenPtsConn
c======================================================================
      implicit none
      include 'vars.h'

      integer          c, i, j, exists, p0, p1, p2, iseed
      double precision x, y, drand, rdx, rdy

c     seed for randum number generator
      iseed = 29374

      do i=1,npmax
         ptype(i)   = -1
         nnbr(i)    = -1
         coord(1,i) = 2.0d0*dmax1(xmin, xmax)
         coord(2,i) = 2.0d0*dmax1(ymin, ymax)
      enddo

      lx = xmax - xmin
      ly = ymax - ymin
      dx = lx/(nx-1)
      dy = ly/(ny-1)
      print*,'lx, ly =',lx,ly
      print*,'dx, dy =',dx,dy

      c  = 0

      print*,'Making coordinates ...'
      do i=1,nx
         do j=1,ny
            c          = c + 1
            x          = xmin + (i-1)*dx
            y          = ymin + (j-1)*dy
            if(perturb.eq.1 .and. i.ne.1 .and. i.ne.nx .and.
     +         j.ne.1 .and. j.ne.ny)then
               rdx = drand(iseed)
               rdy = drand(iseed)
               x   = x + rdx*pertx*dx
               y   = y + rdy*perty*dy
            endif
            coord(1,c) = x
            coord(2,c) = y
            ptype(c)   = INTERIOR
            if(i.eq.nx .or. j.eq.ny)then
               ptype(c)= PERIODIC
            endif
            nnum(i,j) = c
         enddo
      enddo

c     these are points in the computational domain
      np0 = c

c     generate ghost layer
c     bottom layer
      do i=1,nx
         c         = c + 1
         nnum(i,0) = c
         ptype(c)  = PERIODIC
      enddo

c     left layer
      do j=1,ny
         c         = c + 1
         nnum(0,j) = c
         ptype(c)  = PERIODIC
      enddo

c     left-bottom corner point
      c         = c + 1
      nnum(0,0) = c
      ptype(c)  = PERIODIC

      npts = c
      print*,'Number of points =',npts

c     make connectivity
      print*,'Making connectivity ...'
      do i=1,nx-1
         do j=1,ny-1
            c  = 0
            p0 = nnum(i,j)

            if(exists(i-1,j).eq.1)then
               c          = c + 1
               conn(c,p0) = nnum(i-1,j)
            endif

            if(exists(i+1,j).eq.1)then
               c          = c + 1
               conn(c,p0) = nnum(i+1,j)
            endif

            if(exists(i,j-1).eq.1)then
               c          = c + 1
               conn(c,p0) = nnum(i,j-1)
            endif

            if(exists(i,j+1).eq.1)then
               c          = c + 1
               conn(c,p0) = nnum(i,j+1)
            endif

            nnbr(p0) = c

         enddo
      enddo

c     bottom ghost layer
      do i=1,nx-1
         p0         = nnum(i,0)
         p1         = nnum(i,ny-1)
         nnbr(p0)   = 1
         conn(1,p0) = p1
         coord(1,p0)= coord(1,p1)
         coord(2,p0)= coord(2,p1) - ly
      enddo
      p0         = nnum(nx,0)
      p1         = nnum(1,ny-1)
      nnbr(p0)   = 1
      conn(1,p0) = p1
      p2         = nnum(nx,ny-1)
      coord(1,p0)= coord(1,p2)
      coord(2,p0)= coord(2,p2) - ly

c     left ghost layer
      do j=1,ny-1
         p0         = nnum(0,j)
         p1         = nnum(nx-1,j)
         nnbr(p0)   = 1
         conn(1,p0) = p1
         coord(1,p0)= coord(1,p1) - lx
         coord(2,p0)= coord(2,p1)
      enddo
      p0         = nnum(0,ny)
      p1         = nnum(nx-1,1)
      nnbr(p0)   = 1
      conn(1,p0) = p1
      p2         = nnum(nx-1,ny)
      coord(1,p0)= coord(1,p2) - lx
      coord(2,p0)= coord(2,p2)

c     left-botton corner point
      p0         = nnum(0,0)
      p1         = nnum(nx-1,ny-1)
      nnbr(p0)   = 1
      conn(1,p0) = p1
      coord(1,p0)= coord(1,p1) - lx
      coord(2,p0)= coord(2,p1) - ly

c     top periodic layer
      do i=1,nx-1
         p0 = nnum(i,ny)
         nnbr(p0) = 1
         conn(1,p0) = nnum(i,1)
      enddo

c     right periodic layer
      do j=1,ny-1
         p0 = nnum(nx,j)
         nnbr(p0) = 1
         conn(1,p0) = nnum(1,j)
      enddo

c     top-right corner point
      p0 = nnum(nx,ny)
      nnbr(p0) = 1
      conn(1,p0) = nnum(1,1)

c     check
      do i=1,npts
         if(ptype(i).eq.-1)then
            print*,'!!! ptype not set for node =',i
         endif
         if(nnbr(i).eq.-1)then
            print*,'!!! nnbr not set for node =',i
         endif
      enddo

c     make triangles for visualization
      print*,'Making triangles ...'
      c = 0
      do i=1,nx-1
         do j=1,ny-1
            c = c + 1
            tri(1,c) = nnum(i,j)
            tri(2,c) = nnum(i+1,j)
            tri(3,c) = nnum(i+1,j+1)

            c = c + 1
            tri(1,c) = nnum(i,j)
            tri(2,c) = nnum(i+1,j+1)
            tri(3,c) = nnum(i,j+1)
         enddo
      enddo

      ntri = c
      print*,'Number of triangles =',ntri

c     print points for verification
      open(unit=20, file='pts.dat')
      do i=0,nx
         do j=0,ny
            p0 = nnum(i,j)
            write(20,*) coord(1,p0), coord(2,p0)
         enddo
      enddo
      close(20)

      return
      end

c======================================================================
c     point (i,j) exists in the grid
      integer function exists(i, j)
c======================================================================
      implicit none
      include 'vars.h'
      integer i, j

      exists = 1

      if(i.lt.0)  exists = 0
      if(i.gt.nx) exists = 0
      if(j.lt.0)  exists = 0
      if(j.gt.ny) exists = 0

      return
      end

c     -----------------------------------------------------------------
      REAL*8 FUNCTION DRAND(iseed)
c     -----------------------------------------------------------------
c     Selection aleatoire d'un nombre entre 0 et 1 suivant une
c     valeur donnee iseed
c
c                  mod(iseed*7141 + 54773, 259200)
c         ran = -----------------------------------
c                            259200
c     -----------------------------------------------------------------
c
c     Parametres d'appel 
c
      INTEGER iseed
c
c     Variables locales 
c
      INTEGER ia, ic, im
      PARAMETER(ia = 7141, ic = 54773, im = 259200)
c
      iseed    = ABS(MOD(iseed*ia+ic, im))
c
      DRAND    = DBLE(iseed)/DBLE(im)
c
      RETURN
      END
