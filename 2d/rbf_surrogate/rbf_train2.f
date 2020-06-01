      subroutine rbf_train2(nvar, npt, x, y, w, af, rcond)
c----------------------------------------------------------------------------
c*** purpose:
c      determine the weights of a RBF neural network
c*** parameters:
c  nvar             -in-    dimension of design space
c  npt              -in-    number of points
c  x                -in-    training points
c  y                -in-    cost function values
c  w                -out-   table for the metamodel (weights)
c  af               -io-    attenuation factor
c  rcond            -out-   condition number
c
c Note: If you want to optimize af with Rippa method, set af<=0.0 before
c       calling this subroutine.
c----------------------------------------------------------------------------
      implicit none
c*** parameter definition:
      integer          npt,nvar
      double precision x(nvar,*), y(*)
      double precision af, rcond
      double precision w(*)

c*** local variables:
      integer          i, j, ivar, Ipvt(npt)
      double precision dist2, d2(npt,npt), H(npt,npt), Z(npt)
      double precision rippa_af, rippa_af_pso

c     compute the distance between xi and xj; only compute the upper triangular 
c     part and copy to lower part.
      Do i=1,npt
         d2(i,i) = 0.0d0
         Do j=i+1,npt
            dist2 = 0.0d0
            Do ivar=1,nvar
               dist2 = dist2 + ( x(ivar,i) - x(ivar,j) )**2
            EndDo
            d2(i,j) = dist2
         EndDo
         Do j=1,i-1
            d2(i,j) = d2(j,i)
         EndDo
      EndDo

c*** Optimize attenuation factor using Rippa method
      if(af .le. 0.0d0) then
         print*,'Using Rippa method to find attenuation factor'
         af = rippa_af_pso(npt, d2, H, y)
      endif

c*** Building of the matrix for interpolation
      call rbf_matrix(npt, af, d2, H)

c*** Solve the system for interpolation
      call dgeco(H, npt, npt, Ipvt, rcond, Z)
      call dcopy(npt, y, 1, w, 1)
      call dgesl(H, npt, npt, Ipvt, w, 0)

      write(*,'(" Attenuation factor =" f10.4)') af
      write(*,'(" Condition number   =" e14.4)') 1.0d0/rcond
c     print*,'Weights:'
c     write(*,11)(w(i), i=1,npt)
c11    format(e20.10)
      
      return
      end

c----------------------------------------------------------------------------
c Rippa method to optimize attenuation factor
c Minimize Rippa cost function using PSO search
c----------------------------------------------------------------------------
      double precision function rippa_af_pso(npt, d2, H, y)
      implicit none
      include 'Param.h'
      integer          npt
      double precision d2(npt,npt), H(npt,npt), y(npt)

      integer          i, j, iter, MAXITER, npart
      parameter(npart=5)
      double precision af, afmin, afmax, afminl, afmaxl, daf, vmax
      double precision rcond, residue, RESTOL, v(npart), vi, vp, vg,
     +                 vabs, omg, phi1, phi2, hd, alp, vf, afl(npart),
     +                 pafl(npart), cost(npart), pcost(npart), 
     +                 gafl, gcost, r1, r2, dmin, dmax, amin, amax,
     +                 amean
      parameter(MAXITER=100,RESTOL=1.0d-3)
      double precision rippa_cost

      omg  = 1.2
      phi1 = 2.0
      phi2 = 2.0
      hd   = 3
      alp  = 0.99
      vf   = 0.98

C Minimum and maximum distance between particles
      dmin  = 1.0d20
      dmax  =-1.0d20
      do i=1,npt
         do j=i+1,npt
            dmin = dmin1(dmin, d2(i,j))
            dmax = dmax1(dmax, d2(i,j))
         enddo
      enddo
      dmin = dsqrt(dmin)
      dmax = dsqrt(dmax)

      if(dmin.lt.MACHEPS)then
         print*,'===============FATAL==============='
         print*,' Data points are too close'
         print*,'===============FATAL==============='
         stop
      endif

      afmin = dmin
      afmax = dmax
      write(*,10) afmin,afmax
10    format(' Initial range of att. fact. =',2f10.4)

C Initial position of particles
      afminl= dlog10(afmin)
      afmaxl= dlog10(afmax)
      daf   = (afmaxl - afminl)/(npart-1)
      do i=1,npart
         afl(i) = afminl + daf*(i-1)
      enddo

C Initial velocities
      do i=1,npart
         v(i) = 0.0d0
      enddo

C Maximum permissible velocity
      vmax = 0.25d0*(afmaxl - afminl)
      print*,'Maximum permitted velocity =',vmax

C Initial cost and best positions
      gafl  = 0.0d0
      gcost = 1.0d30
      do i=1,npart
         af       = 10.0d0**afl(i)
         cost(i)  = rippa_cost(npt, af, d2, H, y, rcond)
         pcost(i) = cost(i)
         pafl(i)  = afl(i)
         if(cost(i) .lt. gcost)then
            gcost = cost(i)
            gafl  = afl(i)
         endif
      enddo

      print*,'Best initial particle =',10.0d0**gafl
      print*,'   with cost function =',gcost

C PSO iteration loop
      iter    = 0
      residue = RESTOL + 1.0d0
      do while(iter.lt.MAXITER .and. residue.gt.RESTOL)
         amin = 1.0d20
         amax =-1.0d20
         amean= 0.0d0
         do i=1,npart
            r1       = 1.0d0*rand(0)
            r2       = 1.0d0*rand(0)
            vi       = omg*v(i)
            vp       = phi1*r1*(pafl(i) - afl(i))
            vg       = phi2*r2*(gafl    - afl(i))
            v(i)     = 0.729d0*( v(i) + vp + vg )
            vabs     = dabs(v(i))
            if(vabs .gt. vmax)then
               v(i) = (v(i)/vabs)*vmax
            endif
            afl(i)   = afl(i) + v(i)
            af       = 10.0d0**afl(i)
c           if(af .lt. afmin .or. af .gt. afmax) then
            if(af .lt. afmin) then
               afl(i) = afminl + (afmaxl-afminl)*rand(0)
               af     = 10.0d0**afl(i)
               v(i)   = 0.0d0
            endif
            cost(i)  = rippa_cost(npt, af, d2, H, y, rcond)
            if( cost(i) .lt. pcost(i) )then
               pcost(i) = cost(i)
               pafl(i)  = afl(i)
            endif
            amean = amean + af
            amin  = dmin1(amin, af)
            amax  = dmax1(amax, af)
         enddo
         amean = amean/npart

         do i=1,npart
            if(cost(i) .lt. gcost)then
               gcost = cost(i)
               gafl  = afl(i)
            endif
         enddo

         iter = iter + 1
         residue = dabs(amax - amin)/amean
c        write(*,'(i4,3f12.4,2e15.4)') iter, 10.0d0**gafl, amin, amax,
c    +                                gcost, residue
      enddo

      if(iter.ge.MAXITER .and. residue.gt.RESTOL)then
         print*,'======================================================'
         print*,'PSO iterations did not converge to specified tolerance'
         print*,'======================================================'
      endif

C Final value of optimized attenuation factor
      rippa_af_pso = 10.0d0**gafl

      if(rippa_af_pso .lt. afmin .or. rippa_af_pso .gt. afmax) then
         print*,'======================================================'
         print*,'Optimal attenuation factor is outside the min-max rang'
         print*,'======================================================'
      endif

      return
      end

c----------------------------------------------------------------------------
c Rippa cost function
c----------------------------------------------------------------------------
      double precision function rippa_cost(npt, af, d2, H, y, rcond)
      implicit none
      integer          npt
      double precision af, d2(npt,npt), H(npt,npt), y(npt), rcond

      integer          i, j, ipvt(npt)
      double precision err, z(npt), w(npt), det, wi(npt)

      call rbf_matrix(npt, af, d2, H)
      call dgeco(H, npt, npt, Ipvt, rcond, Z)

      if(rcond.lt.1.0d-15)then
         rippa_cost = 1.0d20
      else
         call dcopy(npt, y, 1, w, 1)
         call dgesl(H, npt, npt, Ipvt, w, 0)
c        call dgedi(H, npt, npt, ipvt, det, z, 01)
         rippa_cost = 0.0d0
         do i=1,npt
            do j=1,npt
               wi(j) = 0.0d0
            enddo
            wi(i)=1.0d0
            call dgesl(H, npt, npt, Ipvt, wi, 0)
            err  = w(i)/wi(i)
            rippa_cost = rippa_cost + dabs(err)
         enddo
         rippa_cost = rippa_cost/npt
      endif

      return
      end

c----------------------------------------------------------------------------
c Construct coefficient matrix for RBF network
c Since it is symmetric, construct upper triangular part and copy to lower
c triangular
c  npt              -in-    number of points
c  af               -in-    attenuation factor
c  d2               -in-    distance between data points
c  H                -out-   RBF coeffient matrix
c----------------------------------------------------------------------------
      subroutine rbf_matrix(npt, af, d2, H)
      implicit none
      integer          npt
      double precision af, d2(npt,npt), H(npt,npt)

      integer          i, j
      double precision factor

      factor = 1.0d0/af/af
      do i=1,npt
         H(i,i) = 1.0d0
         do j=i+1,npt
            H(i,j) = dexp(-factor*d2(i,j))
         enddo
         do j=1,i-1
            H(i,j) = H(j,i)
         enddo
      enddo

      return
      end
