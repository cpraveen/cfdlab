c############################################################################
      subroutine rbf_eval2(nvar, npt, x, w, af, xe, fe)
c############################################################################
c*** purpose:
c      determine the weights of a RBF neural network
c*** parameters:
c  npt              -in-    number of points
c  nvar             -in-    dimension of design space
c  x                -in-    training points
c  y                -in-    cost function values
c  w                -out-   table for the metamodel (weights)
c  af               -out-   attenuation factor
c----------------------------------------------------------------------------
      implicit none
c*** parameter definition:
      integer          npt,nvar
      double precision x(nvar,*), xe(nvar), fe
      double precision af
      double precision w(*)

c*** local variables:
      integer          ipt, ivar
      double precision factor, dist2

      fe     = 0.0d0
      factor = 1.0d0/af/af
      do ipt=1,npt
         dist2 = 0.0d0
         do ivar=1,nvar
            dist2 = dist2 + ( xe(ivar) - x(ivar,ipt) )**2
         enddo
         fe = fe + w(ipt)*dexp(-factor*dist2)
      enddo

      return
      end

c############################################################################
      subroutine rbf_eval_grad(nvar, npt, x, w, af, xe, fed)
c############################################################################
c*** purpose:
c      determine the weights of a RBF neural network
c*** parameters:
c  npt              -in-    number of points
c  nvar             -in-    dimension of design space
c  x                -in-    training points
c  y                -in-    cost function values
c  w                -out-   table for the metamodel (weights)
c  af               -out-   attenuation factor
c----------------------------------------------------------------------------
      implicit none
c*** parameter definition:
      integer          npt,nvar
      double precision x(nvar,*), xe(nvar), fed(nvar)
      double precision af
      double precision w(*)

c*** local variables:
      integer          ipt, ivar
      double precision factor, fact2, dist2, dx(nvar)

      do ivar=1,nvar
         fed(ivar)     = 0.0d0
      enddo

      factor = 1.0d0/af/af
      do ipt=1,npt
         dist2 = 0.0d0
         do ivar=1,nvar
            dx(ivar) = xe(ivar) - x(ivar,ipt)
            dist2 = dist2 + dx(ivar)**2
         enddo
         fact2 = -2.0d0*factor*w(ipt)*dexp(-factor*dist2)
         do ivar=1,nvar
            fed(ivar) = fed(ivar) + fact2*dx(ivar)
         enddo
      enddo

      return
      end

c############################################################################
      subroutine rbf_eval_hess(nvar, npt, x, w, af, xe, fedd)
c############################################################################
c*** purpose:
c      determine hessian at xe using RBF neural network
c*** parameters:
c  npt              -in-    number of points
c  nvar             -in-    dimension of design space
c  x                -in-    training points
c  y                -in-    cost function values
c  w                -out-   table for the metamodel (weights)
c  af               -out-   attenuation factor
c----------------------------------------------------------------------------
      implicit none
c*** parameter definition:
      integer          npt,nvar
      double precision x(nvar,*), xe(nvar), fedd(nvar,nvar)
      double precision af
      double precision w(*)

c*** local variables:
      integer          ipt, ivar, jvar
      double precision factor, fact2, dist2, dx(nvar)

c***  set zero on upper triangular part including on diagonal
      do ivar=1,nvar
         do jvar=ivar,nvar
            fedd(ivar,jvar) = 0.0d0
         enddo
      enddo

      factor = 1.0d0/af/af
      do ipt=1,npt
         dist2 = 0.0d0
         do ivar=1,nvar
            dx(ivar) = xe(ivar) - x(ivar,ipt)
            dist2 = dist2 + dx(ivar)**2
         enddo
         fact2 = -2.0d0*factor*w(ipt)*dexp(-factor*dist2)
         do ivar=1,nvar
            fedd(ivar,ivar) = fedd(ivar,ivar) + fact2*(1.0d0 - 
     +                        2.0d0*factor*dx(ivar)**2)
            do jvar=ivar+1,nvar
               fedd(ivar,jvar) = fedd(ivar,jvar) -
     +                           fact2*factor*dx(ivar)*dx(jvar)
            enddo
         enddo
c        copy upper triangular part to lower triangular
         do ivar=2,nvar
            do jvar=1,ivar-1
               fedd(ivar,jvar) = fedd(jvar,ivar)
            enddo
         enddo
      enddo

      return
      end
