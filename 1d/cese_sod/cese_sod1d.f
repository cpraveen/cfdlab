c                      Sod's 1-D shock tube problem
c                        A Sample Program from:
c
c     S.C. Chang, S.T. Yu, A. Himansu, X.Y. Wang, C.Y. Chow and C.Y. Loh, 
c     "The Method of Space-Time Conservation Element and Solution Element -- 
c     A New Paradigm for Numerical Solution of Conservation Laws," 
c     Computational Fluid Dynamics Review 1998, M. Hafez and K. Oshima, eds., 
c     World Scientific, New Jersey, USA. 

      implicit real*8(a-h,o-z)
      parameter (nxd=1000)
      dimension  q(3,nxd), qn(3,nxd), qx(3,nxd), qt(3,nxd),
     .           s(3,nxd), vxl(3), vxr(3), xx(nxd)

c     nx must be an odd integer.
      nx = 101
      it = 100
      dt = 0.4d-2
      dx = 0.1d-1
      ga = 1.4d0
      rhol = 1.d0
      ul = 0.d0
      pl = 1.d0
      rhor = 0.125d0
      ur = 0.d0
      pr = 0.1d0
      ia = 2

      nx1 = nx + 1
      nx2 = nx1/2
      hdt = dt/2.d0
      tt = hdt*dfloat(it)
      qdt = dt/4.d0
      hdx = dx/2.d0
      qdx = dx/4.d0
      dtx = dt/dx
      a1 = ga - 1.d0
      a2 = 3.d0 - ga
      a3 = a2/2.d0
      a4 = 1.5d0*a1
      u2l = rhol*ul
      u3l = pl/a1 + 0.5d0*rhol*ul**2
      u2r = rhor*ur
      u3r = pr/a1 + 0.5d0*rhor*ur**2

      do j = 1,nx2
       q(1,j) = rhol
       q(2,j) = u2l
       q(3,j) = u3l
       q(1,nx2+j) = rhor
       q(2,nx2+j) = u2r
       q(3,nx2+j) = u3r
       do i = 1,3
        qx(i,j) = 0.d0
        qx(i,nx2+j) = 0.d0
       enddo
      enddo

      open (unit=8,file='for008')
      write (8,10) tt,it,ia,nx
      write (8,20) dt,dx,ga
      write (8,30) rhol,ul,pl
      write (8,40) rhor,ur,pr

      do i = 1,it
       m = nx + i - (i/2)*2
       do j = 1,m
        w2 = q(2,j)/q(1,j)
        w3 = q(3,j)/q(1,j)
        f21 = -a3*w2**2
        f22 = a2*w2
        f31 = a1*w2**3 - ga*w2*w3
        f32 = ga*w3 - a4*w2**2
        f33 = ga*w2
        qt(1,j) = -qx(2,j)
        qt(2,j) = -(f21*qx(1,j) + f22*qx(2,j) + a1*qx(3,j))
        qt(3,j) = -(f31*qx(1,j) + f32*qx(2,j) + f33*qx(3,j))
        s(1,j) = qdx*qx(1,j) + dtx*(q(2,j) + qdt*qt(2,j))
        s(2,j) = qdx*qx(2,j) + dtx*(f21*(q(1,j) + qdt*qt(1,j)) +
     .    f22*(q(2,j) + qdt*qt(2,j)) + a1*(q(3,j) + qdt*qt(3,j)))
        s(3,j) = qdx*qx(3,j) + dtx*(f31*(q(1,j) + qdt*qt(1,j)) +
     .    f32*(q(2,j) + qdt*qt(2,j)) + f33*(q(3,j) + qdt*qt(3,j)))
       enddo

       if (i.ne.(i/2)*2) goto 150
        do k = 1,3
         qx(k,nx1) = qx(k,nx)
         qn(k,1) = q(k,1)
         qn(k,nx1) = q(k,nx)
        enddo

150    j1 = 1 - i + (i/2)*2
       mm = m - 1
       do j = 1,mm
        do k = 1,3
         qn(k,j+j1) = 0.5d0*(q(k,j) + q(k,j+1) + s(k,j) - s(k,j+1))
         vxl(k) = (qn(k,j+j1) - q(k,j) - hdt*qt(k,j))/hdx
         vxr(k) = (q(k,j+1) + hdt*qt(k,j+1) - qn(k,j+j1))/hdx
         qx(k,j+j1) = (vxl(k)*(dabs(vxr(k)))**ia + vxr(k)*(dabs(vxl(k)))
     .       **ia)/((dabs(vxl(k)))**ia + (dabs(vxr(k)))**ia + 1.d-60)
        enddo
       enddo

      m = nx1 - i + (i/2)*2
       do j = 1,m
        do k = 1,3
         q(k,j) = qn(k,j)
        enddo
       enddo
      enddo

      m = nx1 -it + (it/2)*2
      mm = m - 1
      xx(1) = -0.5d0*dx*dfloat(mm)
      do j = 1,mm
       xx(j+1) = xx(j) + dx
      enddo
      do j = 1,m
       x = q(2,j)/q(1,j)
       y = a1*(q(3,j) - 0.5d0*x**2*q(1,j))
       z = x/dsqrt(ga*y/q(1,j))
       write (8,50) xx(j),q(1,j),x,y,z
       write (9,55) xx(j),q(1,j),x,y,z
      enddo

      close (unit=8)
10    format(' t = ',g14.7,' it = ',i4,' ia = ',i4,' nx = ',i4)
20    format(' dt = ',g14.7,' dx = ',g14.7,' gamma = ',g14.7)
30    format(' rhol = ',g14.7,' ul = ',g14.7,' pl = ',g14.7)
40    format(' rhor = ',g14.7,' ur = ',g14.7,' pr = ',g14.7)
50    format(' x =',f8.4,' rho =',f8.4,' u =',f8.4,' p =',f8.4,
     .       ' M =',f8.4)
55    format(5f8.4)

      stop
      end