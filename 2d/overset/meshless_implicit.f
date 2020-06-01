c======================================================================
      subroutine meshless_implicit(dt, Q0, Q1, res, cx, cy, dQ0, dQ1)
c======================================================================
c     output: dQ0
c     input : everything else
c     NOTE: Q0, Q1 are primitive variables but dQ0, dQ1 are in terms of
c     conserved variables
c     dt             : time step
c     Q0(nvar)       : primitive vars
c     Q1(nvar,nnbr)  : primitive vars
c     res(nvar)      : meshless residual computed in meshless.f
c     cx, cy         : coefficients in least-squares formula
c     dQ0(nvar)      : change in conserved vars
c     dQ1(nvar,nnbr) : change in conserved vars
c======================================================================
      implicit none
      include 'misc.h'
      real    dt, Q0(*), Q1(nvar,*), res(*), cx(*), cy(*), dQ0(*), 
     +        dQ1(nvar,*)

      integer i, j, k, ipvt(nvar), info
      real    jac0(nvar,nvar), jac1(nvar,nvar), jac2(nvar,nvar),
     +        dmat(nvar,nvar), Alhs(nvar,nvar), rhs(nvar)

      Alhs(1:nvar,1:nvar) = 0.0
      rhs(1:nvar) = -res(1:nvar)

      do i=1,nnbr
         call euler_jacobian(Q0,      cx(i), cy(i), jac0)
         call euler_jacobian(Q1(1,i), cx(i), cy(i), jac1)
         call diss_matrix(Q0, Q1(1,i),cx(i), cy(i), dmat)

c        right hand side
         jac2(1:nvar,1:nvar) = jac1(1:nvar,1:nvar) - dmat(1:nvar,1:nvar)
         do j=1,nvar
            do k=1,nvar
               rhs(j) = rhs(j) + jac2(j,k)*dQ1(k,i)
            enddo
         enddo

c        left hand side
         Alhs(1:nvar,1:nvar) = Alhs(1:nvar,1:nvar) - jac0(1:nvar,1:nvar)
     +                               + dmat(1:nvar,1:nvar)
      enddo

      do i=1,nvar
         Alhs(i,i) = 1.0/dt + Alhs(i,i)
      enddo

c     solve using linpack Alhs * dQ0 = rhs
      call sgefa(Alhs, nvar, nvar, ipvt, info)
      if(info.ne.0)then
         print*,'meshless_implicit: error, info is non-zero'
         print*,'info =',info
         stop
      endif
      call sgesl(Alhs, nvar, nvar, ipvt, rhs, 0)
      dQ0(1:nvar) = rhs(1:nvar)

      return
      end

c======================================================================
      subroutine euler_jacobian(q, nx, ny, jac)
c======================================================================
      implicit none
      real q(*), nx, ny, jac(4,4)

      real gm, gm1, gm3
      real r, u, v, p, un, c2, h, u2

      gm = 1.4
      gm1= gm - 1.0
      gm3= 3.0 - gm

      r = q(1)
      u = q(2)
      v = q(3)
      p = q(4)

      un= u*nx + v*ny
      u2= u**2 + v**2
      c2= gm*p/r
      h = c2/gm1 + 0.5*u2

      jac(1,1) = 0.0
      jac(2,1) = -u*un + 0.5*gm1*u2*nx
      jac(3,1) = -v*un + 0.5*gm1*u2*ny
      jac(4,1) = -(h - 0.5*gm1*u2)*un

      jac(1,2) = nx
      jac(2,2) = gm3*u*nx + v*ny
      jac(3,2) = -gm1*u*ny + v*nx
      jac(4,2) = h*nx - gm1*u*un

      jac(1,3) = ny
      jac(2,3) = -gm1*v*nx + u*ny
      jac(3,3) = gm3*v*ny + u*nx
      jac(4,3) = h*ny - gm1*v*un

      jac(1,4) = 0.0
      jac(2,4) = gm1*nx
      jac(3,4) = gm1*ny
      jac(4,4) = gm*un

      return
      end

C------------------------------------------------------------------------------
C     Computes euler split flux jacobians
C------------------------------------------------------------------------------
      subroutine diss_matrix(q1, q2, nx, ny, jac)
      implicit none
      real    q1(*), q2(*), nx, ny, jac(4,4)

      integer i, j, k
      real    r1,u1,v1,p1,h1,c12
      real    r2,u2,v2,p2,h2,c22
      real    sr1,sr2,srd,vel2
      real    gm, gm1
      real    ua, va, pa, qa2, aa2, aa, ha, ra,
     &        una, vna, ct, st, lent,
     &        m2, t1, t2, t3, t4, t5, l1, l2, l3, l4,
     &        S(4,4), R(4,4), jtmp

c     gm1 = gamma - 1
      gm = 1.4
      gm1= gm - 1.0

      lent = sqrt(nx**2 + ny**2)
      ct = nx/lent
      st = ny/lent

c     left state
      r1 = q1(1)
      u1 = q1(2)
      v1 = q1(3)
      p1 = q1(4)
      c12= gm*p1/r1
      vel2 = u1**2 + v1**2
      h1 = c12/gm1 + 0.5*vel2

c     right state
      r2 = q2(1)
      u2 = q2(2)
      v2 = q2(3)
      p2 = q2(4)
      c22= gm*p2/r2
      vel2 = u2**2 + v2**2
      h2 = c22/gm1 + 0.5*vel2

C     roe-average State
      sr1= sqrt(r1)
      sr2= sqrt(r2)
      srd= 1.0/(sr1 + sr2)
      ra = sqrt(r1*r2)
      ua = (sr1*u1 + sr2*u2)*srd
      va = (sr1*v1 + sr2*v2)*srd
      ha = (sr1*h1 + sr2*h2)*srd
      qa2= ua**2 + va**2
      aa2= gm1*(ha - 0.5*qa2)
      aa = sqrt(aa2)
      pa = aa2*ra/gm

C     Rotated velocity
      una = ua*ct + va*st
      vna =-ua*st + va*ct

C     absolute Eigenvalues
      l1 = abs(una)
      l2 = l1
      l3 = abs(una + aa)
      l4 = abs(una - aa)

c     Right eigenvector matrix
      t1 = 0.5*ra/aa
      m2 = qa2/aa2

      R(1,1) = 1.0
      R(2,1) = ua
      R(3,1) = va
      R(4,1) = 0.5*qa2

      R(1,2) = 0.0
      R(2,2) = ra*st
      R(3,2) = -ra*ct
      R(4,2) = -ra*vna

      R(1,3) = t1
      R(2,3) = t1*(ua + aa*ct)
      R(3,3) = t1*(va + aa*st)
      R(4,3) = t1*(ha + aa*una)

      R(1,4) = t1
      R(2,4) = t1*(ua - aa*ct)
      R(3,4) = t1*(va - aa*st)
      R(4,4) = t1*(ha - aa*una)

c     Inverse of right eigenvector matrix
      t1     = 0.5*gm1*m2*aa/ra
      t2     = una/ra
      t3     = gm1*ua/aa
      t4     = gm1*va/aa
      t5     = gm1/(ra*aa)

      S(1,1) = 1.0 - 0.5*gm1*m2
      S(2,1) = vna/ra
      S(3,1) = t1 - t2
      S(4,1) = t1 + t2

      S(1,2) = t3/aa
      S(2,2) = st/ra
      S(3,2) = (ct - t3)/ra
      S(4,2) =-(ct + t3)/ra

      S(1,3) = t4/aa
      S(2,3) = -ct/ra
      S(3,3) = (st - t4)/ra
      S(4,3) =-(st + t4)/ra

      S(1,4) = -gm1/aa2
      S(2,4) = 0.0
      S(3,4) = t5
      S(4,4) = t5

c     Multiply R * lambda
      do i=1,4
         R(i,1) = l1*R(i,1)
         R(i,2) = l2*R(i,2)
         R(i,3) = l3*R(i,3)
         R(i,4) = l4*R(i,4)
      enddo

c     Now multiply with S
      do i=1,4
         do j=1,4
            jtmp = 0.0
            do k=1,4
               jtmp = jtmp + R(i,k)*S(k,j)
            enddo
            jac(i,j) = lent*jtmp
         enddo
      enddo

      return
      end

c======================================================================
c routines from linpack for solution of Ax=b
c======================================================================
      subroutine sgefa(a,lda,n,ipvt,info)
c======================================================================
      integer lda,n,ipvt(1),info
      real a(lda,1)
c
c     sgefa factors a real matrix by gaussian elimination.
c
c     sgefa is usually called by sgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for sgeco) = (1 + 9/n)*(time for sgefa) .
c
c     on entry
c
c        a       real(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that sgesl or sgedi will divide by zero
c                     if called.  use  rcond  in sgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas saxpy,sscal,isamax
c
c     internal variables
c
      real t
      integer isamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = isamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0e0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0e0/a(k,k)
            call sscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call saxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0e0) info = n
      return
      end

c======================================================================
      subroutine sgesl(a,lda,n,ipvt,b,job)
c======================================================================
      integer lda,n,ipvt(1),job
      real a(lda,1),b(1)
c
c     sgesl solves the real system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by sgeco or sgefa.
c
c     on entry
c
c        a       real(lda, n)
c                the output from sgeco or sgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from sgeco or sgefa.
c
c        b       real(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if sgeco has set rcond .gt. 0.0
c        or sgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call sgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call sgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas saxpy,sdot
c
c     internal variables
c
      real sdot,t
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call saxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call saxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            t = sdot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + sdot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end

c======================================================================
      SUBROUTINE SSCAL(N,SA,SX,INCX)
c======================================================================
*     .. Scalar Arguments ..
      REAL SA
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      REAL SX(*)
*     ..
*
*  Purpose
*  =======
*
*     scales a vector by a constant.
*     uses unrolled loops for increment equal to 1.
*     jack dongarra, linpack, 3/11/78.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
      INTEGER I,M,MP1,NINCX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 20
*
*        code for increment not equal to 1
*
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
          SX(I) = SA*SX(I)
   10 CONTINUE
      RETURN
*
*        code for increment equal to 1
*
*
*        clean-up loop
*
   20 M = MOD(N,5)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          SX(I) = SA*SX(I)
   30 CONTINUE
      IF (N.LT.5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
          SX(I) = SA*SX(I)
          SX(I+1) = SA*SX(I+1)
          SX(I+2) = SA*SX(I+2)
          SX(I+3) = SA*SX(I+3)
          SX(I+4) = SA*SX(I+4)
   50 CONTINUE
      RETURN
      END

c======================================================================
      SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)
c======================================================================
*     .. Scalar Arguments ..
      REAL SA
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      REAL SX(*),SY(*)
*     ..
*
*  Purpose
*  =======
*
*     SAXPY constant times a vector plus a vector.
*     uses unrolled loop for increments equal to one.
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0) RETURN
      IF (SA.EQ.0.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
*
*        code for unequal increments or equal increments
*          not equal to 1
*
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          SY(IY) = SY(IY) + SA*SX(IX)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
   20 M = MOD(N,4)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          SY(I) = SY(I) + SA*SX(I)
   30 CONTINUE
      IF (N.LT.4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
          SY(I) = SY(I) + SA*SX(I)
          SY(I+1) = SY(I+1) + SA*SX(I+1)
          SY(I+2) = SY(I+2) + SA*SX(I+2)
          SY(I+3) = SY(I+3) + SA*SX(I+3)
   50 CONTINUE
      RETURN
      END

c======================================================================
      REAL FUNCTION SDOT(N,SX,INCX,SY,INCY)
c======================================================================
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      REAL SX(*),SY(*)
*     ..
*
*  Purpose
*  =======
*
*     forms the dot product of two vectors.
*     uses unrolled loops for increments equal to one.
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*

*     .. Local Scalars ..
      REAL STEMP
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      STEMP = 0.0e0
      SDOT = 0.0e0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
*
*        code for unequal increments or equal increments
*          not equal to 1
*
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          STEMP = STEMP + SX(IX)*SY(IY)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      SDOT = STEMP
      RETURN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
   20 M = MOD(N,5)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          STEMP = STEMP + SX(I)*SY(I)
   30 CONTINUE
      IF (N.LT.5) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
          STEMP = STEMP + SX(I)*SY(I) + SX(I+1)*SY(I+1) +
     +            SX(I+2)*SY(I+2) + SX(I+3)*SY(I+3) + SX(I+4)*SY(I+4)
   50 CONTINUE
   60 SDOT = STEMP
      RETURN
      END

c======================================================================
      INTEGER FUNCTION ISAMAX(N,SX,INCX)
c======================================================================
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      REAL SX(*)
*     ..
*
*  Purpose
*  =======
*
*     finds the index of element having max. absolute value.
*     jack dongarra, linpack, 3/11/78.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
      REAL SMAX
      INTEGER I,IX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC ABS
*     ..
      ISAMAX = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      ISAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) GO TO 20
*
*        code for increment not equal to 1
*
      IX = 1
      SMAX = ABS(SX(1))
      IX = IX + INCX
      DO 10 I = 2,N
          IF (ABS(SX(IX)).LE.SMAX) GO TO 5
          ISAMAX = I
          SMAX = ABS(SX(IX))
    5     IX = IX + INCX
   10 CONTINUE
      RETURN
*
*        code for increment equal to 1
*
   20 SMAX = ABS(SX(1))
      DO 30 I = 2,N
          IF (ABS(SX(I)).LE.SMAX) GO TO 30
          ISAMAX = I
          SMAX = ABS(SX(I))
   30 CONTINUE
      RETURN
      END
