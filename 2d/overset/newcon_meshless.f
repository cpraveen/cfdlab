c----------------------------------------------------------------------------
      subroutine  do_interpolations_meshless(q,vnu,x,y,tscale,jmax,kmax,
     <     q1,vnu1,x1,y1,xx1,xy1,yx1,yy1,jmax1,kmax1,
     <     imesh,idonor,frac,ndonor,Mpts)
c----------------------------------------------------------------------------
c---- Use meshless method to update chimera points.
c---- q,vnu,x,y      : Mesh that contains chimera  points
c---- q1,vnu1,x1,y1  : Mesh that contributes donor points
c---------------------------------------------------------------------------

      implicit none

      integer jmax,kmax,jmax1,kmax1
      real q(jmax,kmax,5),x(jmax,kmax),y(jmax,kmax),tscale(jmax,kmax)
      real q1(jmax1,kmax1,5),x1(jmax1,kmax1),y1(jmax1,kmax1)
      real xx1(jmax1,kmax1),xy1(jmax1,kmax1)
      real yx1(jmax1,kmax1),yy1(jmax1,kmax1)
      real vnu(jmax,kmax),vnu1(jmax1,kmax1)
      integer imesh(Mpts,3),idonor(Mpts,3)
      real frac(Mpts,3),frac1(Mpts,3)
      real gm1,rhoi,dj,dk

      real, allocatable :: x_0(:),y_0(:),x_1(:,:),y_1(:,:),dt(:)
      real, allocatable :: Q_0(:,:),Q_1(:,:,:),Qx_1(:,:,:),Qy_1(:,:,:)

      integer id,ndonor,Mpts,i,j,nbr,nnbr
      integer ii,jj,iipp,jjpp,iip,jjp,iim,jjm,i1,j1,ipoint(4),jpoint(4)
      integer ip1,im1,jp1,jm1
      real qip1,qip2,qip3,qip4,qim1,qim2,qim3,qim4
      real qjp1,qjp2,qjp3,qjp4,qjm1,qjm2,qjm3,qjm4
 
      gm1=0.4
      nnbr=4

      allocate(x_0(ndonor),y_0(ndonor),dt(ndonor),Q_0(4,ndonor))
      allocate(x_1(nnbr,ndonor),y_1(nnbr,ndonor),Q_1(4,nnbr,ndonor))
      allocate(Qx_1(4,nnbr,ndonor),Qy_1(4,nnbr,ndonor))


      do id=1,ndonor

c	 Get pointers --> (i,j) is point in question, 
c                         (ii,jj) is point in donor mesh
         
         i=imesh(id,1)
         j=imesh(id,2)
         ii=idonor(id,1)
         jj=idonor(id,2)

         iipp=min(ii+2,jmax1)
         jjpp=min(jj+2,kmax1)
         iip=min(ii+1,jmax1)
         jjp=min(jj+1,kmax1)
         iim=max(ii-1,1)
         jjm=min(jj-1,1)

c         Populate array of "Unknowns"
 
         dt(id) = tscale(i,j)

         x_0(id) = x(i,j)
         y_0(id) = y(i,j)

          rhoi      = 1.0/q(i,j,1)
          Q_0(1,id) = q(i,j,1)*q(i,j,5)
          Q_0(2,id) = q(i,j,2)*rhoi
          Q_0(3,id) = q(i,j,3)*rhoi
          Q_0(4,id) = gm1*( q(i,j,4)-0.5*(q(i,j,2)**2+q(i,j,3)**2)
     <                                              /q(i,j,1))*q(i,j,5)

c         Get pointers to fill array of "Knowns"

          ipoint(1)=ii
          ipoint(2)=iip
          ipoint(3)=iip
          ipoint(4)=ii

          jpoint(1)=jj
          jpoint(2)=jj
          jpoint(3)=jjp
          jpoint(4)=jjp

c         Ok.. those were under circumstances where the neighbors were decently distributed
c         If not, widen stencil

          dj=frac(id,1)
          dk=frac(id,2)

          if(dj.lt.0.25) then
          ipoint(1)=iim
          ipoint(4)=iim
          endif

          if(dj.gt.0.75) then
          ipoint(2)=iipp
          ipoint(3)=iipp
          endif

          if(dk.lt.0.25) then
          jpoint(1)=jjm
          jpoint(2)=jjm
          endif

          if(dk.gt.0.75) then
          jpoint(3)=jjpp
          jpoint(4)=jjpp
          endif


c         Now populate array of "Knowns"

          do nbr=1,4
 
          i1=ipoint(nbr)
          j1=jpoint(nbr)
          ip1=min(i1+1,jmax1)
          im1=max(i1-1,1)
          jp1=min(j1+1,kmax1)
          jm1=max(j1-1,1)

          x_1(nbr,id)  = x1(i1,j1)
          y_1(nbr,id)  = y1(i1,j1)
          rhoi         = 1.0/q1(i1,j1,1)
          Q_1(1,nbr,id) = q1(i1,j1,1)*q1(i1,j1,5)
          Q_1(2,nbr,id) = q1(i1,j1,2)*rhoi
          Q_1(3,nbr,id) = q1(i1,j1,3)*rhoi
          Q_1(4,nbr,id) = gm1*( q1(i1,j1,4)-0.5*(q1(i1,j1,2)**2+q1(i1,j1,3)**2)
     <                                         /q1(i1,j1,1))*q1(i1,j1,5)

c         Gradients

          rhoi          = 1.0/q1(ip1,j1,1)
          qip1 		= q1(ip1,j1,1)*q1(ip1,j1,5)
          qip2          = q1(ip1,j1,2)*rhoi
          qip3		= q1(ip1,j1,3)*rhoi
          qip4          = gm1*( q1(ip1,j1,4)-0.5*(q1(ip1,j1,2)**2+q1(ip1,j1,3)**2)
     <                                         /q1(ip1,j1,1))*q1(ip1,j1,5)
          rhoi          = 1.0/q1(im1,j1,1)
          qim1 		= q1(im1,j1,1)*q1(im1,j1,5)
          qim2          = q1(im1,j1,2)*rhoi
          qim3		= q1(im1,j1,3)*rhoi
          qim4          = gm1*( q1(im1,j1,4)-0.5*(q1(im1,j1,2)**2+q1(im1,j1,3)**2)
     <                                         /q1(im1,j1,1))*q1(im1,j1,5)

          rhoi          = 1.0/q1(i1,jp1,1)
          qjp1 		= q1(i1,jp1,1)*q1(i1,jp1,5)
          qjp2          = q1(i1,jp1,2)*rhoi
          qjp3		= q1(i1,jp1,3)*rhoi
          qjp4          = gm1*( q1(i1,jp1,4)-0.5*(q1(i1,jp1,2)**2+q1(i1,jp1,3)**2)
     <                                         /q1(i1,jp1,1))*q1(i1,jp1,5)
          rhoi          = 1.0/q1(i1,jm1,1)
          qjm1 		= q1(i1,jm1,1)*q1(i1,jm1,5)
          qjm2          = q1(i1,jm1,2)*rhoi
          qjm3		= q1(i1,jm1,3)*rhoi
          qjm4          = gm1*( q1(i1,jm1,4)-0.5*(q1(i1,jm1,2)**2+q1(i1,jm1,3)**2)
     <                                         /q1(i1,jm1,1))*q1(i1,jm1,5)

          Qx_1(1,nbr,id)=0.5*((qip1-qim1)*xx1(i1,j1)+(qjp1-qjm1)*yx1(i1,j1))
          Qx_1(2,nbr,id)=0.5*((qip2-qim2)*xx1(i1,j1)+(qjp2-qjm2)*yx1(i1,j1))
          Qx_1(3,nbr,id)=0.5*((qip3-qim3)*xx1(i1,j1)+(qjp3-qjm3)*yx1(i1,j1))
          Qx_1(4,nbr,id)=0.5*((qip4-qim4)*xx1(i1,j1)+(qjp4-qjm4)*yx1(i1,j1))

          Qy_1(1,nbr,id)=0.5*((qip1-qim1)*xy1(i1,j1)+(qjp1-qjm1)*yy1(i1,j1))
          Qy_1(2,nbr,id)=0.5*((qip2-qim2)*xy1(i1,j1)+(qjp2-qjm2)*yy1(i1,j1))
          Qy_1(3,nbr,id)=0.5*((qip3-qim3)*xy1(i1,j1)+(qjp3-qjm3)*yy1(i1,j1))
          Qy_1(4,nbr,id)=0.5*((qip4-qim4)*xy1(i1,j1)+(qjp4-qjm4)*yy1(i1,j1))


          enddo

      enddo

c        Send stuff to Praveendora's Magic Box

         call meshless(dt,ndonor,x_0,y_0,Q_0,x_1,y_1,Q_1,Qx_1,Qy_1)

c        Update actually happens here. Q_0 is just dt*Residue

      do id=1,ndonor
         
         i=imesh(id,1)
         j=imesh(id,2)

         q(i,j,1)=q(i,j,1)+Q_0(1,id)/q(i,j,5)
         q(i,j,2)=q(i,j,2)+Q_0(2,id)/q(i,j,5)
         q(i,j,3)=q(i,j,3)+Q_0(3,id)/q(i,j,5)
         q(i,j,4)=q(i,j,4)+Q_0(4,id)/q(i,j,5)

      enddo

      return
      end


c----------------------------------------------------------------------------
      subroutine  do_interpolations_meshless_high(q,vnu,x,y,tscale,jmax,kmax,
     <     q1,vnu1,x1,y1,xx1,xy1,yx1,yy1,jmax1,kmax1,
     <     imesh,idonor,frac,ndonor,Mpts)
c----------------------------------------------------------------------------
c---- Use meshless method to update chimera points.
c---- q,vnu,x,y      : Mesh that contains chimera  points
c---- q1,vnu1,x1,y1  : Mesh that contributes donor points
c---------------------------------------------------------------------------

      implicit none

      integer jmax,kmax,jmax1,kmax1
      real q(jmax,kmax,5),x(jmax,kmax),y(jmax,kmax),tscale(jmax,kmax)
      real q1(jmax1,kmax1,5),x1(jmax1,kmax1),y1(jmax1,kmax1)
      real xx1(jmax1,kmax1),xy1(jmax1,kmax1)
      real yx1(jmax1,kmax1),yy1(jmax1,kmax1)
      real vnu(jmax,kmax),vnu1(jmax1,kmax1)
      integer imesh(Mpts,3),idonor(Mpts,3)
      real frac(Mpts,3),frac1(Mpts,3)
      real gm1,rhoi,dj,dk

      real, allocatable :: x_0(:),y_0(:),x_1(:,:),y_1(:,:),dt(:)
      real, allocatable :: Q_0(:,:),Q_1(:,:,:),Qx_1(:,:,:),Qy_1(:,:,:)

      integer, allocatable :: ipoint(:),jpoint(:)

      integer id,ndonor,Mpts,i,j,nbr,nnbr
      integer ii,jj,iipp,jjpp,iip,jjp,iim,jjm,i1,j1
      integer ip1,im1,jp1,jm1
      real qip1,qip2,qip3,qip4,qim1,qim2,qim3,qim4
      real qjp1,qjp2,qjp3,qjp4,qjm1,qjm2,qjm3,qjm4
 
      gm1=0.4
      nnbr=12

      allocate(x_0(ndonor),y_0(ndonor),dt(ndonor),Q_0(4,ndonor))
      allocate(x_1(nnbr,ndonor),y_1(nnbr,ndonor),Q_1(4,nnbr,ndonor))
      allocate(Qx_1(4,nnbr,ndonor),Qy_1(4,nnbr,ndonor))
      allocate(ipoint(nnbr),jpoint(nnbr))


      do id=1,ndonor

c	 Get pointers --> (i,j) is point in question, 
c                         (ii,jj) is point in donor mesh
         
         i=imesh(id,1)
         j=imesh(id,2)
         ii=idonor(id,1)
         jj=idonor(id,2)

         iipp=min(ii+2,jmax1)
         jjpp=min(jj+2,kmax1)
         iip=min(ii+1,jmax1)
         jjp=min(jj+1,kmax1)
         iim=max(ii-1,1)
         jjm=min(jj-1,1)

c        Populate array of "Unknowns"
 
         dt(id) = tscale(i,j)

         x_0(id) = x(i,j)
         y_0(id) = y(i,j)

          rhoi      = 1.0/q(i,j,1)
          Q_0(1,id) = q(i,j,1)*q(i,j,5)
          Q_0(2,id) = q(i,j,2)*rhoi
          Q_0(3,id) = q(i,j,3)*rhoi
          Q_0(4,id) = gm1*( q(i,j,4)-0.5*(q(i,j,2)**2+q(i,j,3)**2)
     <                                              /q(i,j,1))*q(i,j,5)

c         Get pointers to fill array of "Knowns"

          ipoint(1)=ii
          ipoint(2)=iip
          ipoint(3)=iip
          ipoint(4)=ii

          jpoint(1)=jj
          jpoint(2)=jj
          jpoint(3)=jjp
          jpoint(4)=jjp

c         Ok.. those were under circumstances where the neighbors were decently distributed
c         If not, widen stencil

          dj=frac(id,1)
          dk=frac(id,2)

          if(dj.lt.0.25) then
          ipoint(1)=iim
          ipoint(4)=iim
          endif

          if(dj.gt.0.75) then
          ipoint(2)=iipp
          ipoint(3)=iipp
          endif

          if(dk.lt.0.25) then
          jpoint(1)=jjm
          jpoint(2)=jjm
          endif

          if(dk.gt.0.75) then
          jpoint(3)=jjpp
          jpoint(4)=jjpp
          endif

          ipoint(5) =ipoint(1)
          ipoint(6) =ipoint(2)
          ipoint(7) =min(ipoint(2)+1,jmax1)
          ipoint(8) =ipoint(7)
          ipoint(9) =ipoint(3)
          ipoint(10)=ipoint(4)
          ipoint(11)=max(ipoint(4)-1,1)
          ipoint(12)=ipoint(11)

          jpoint(5) =max(jpoint(1)-1,1)
          jpoint(6) =jpoint(5)
          jpoint(7) =jpoint(2)
          jpoint(8) =jpoint(3)
          jpoint(9) =min(jpoint(3)+1,kmax1)
          jpoint(10)=jpoint(9)
          jpoint(11)=jpoint(4)
          jpoint(12)=jpoint(1)


c         Now populate array of "Knowns"

          do nbr=1,nnbr
 
          i1=ipoint(nbr)
          j1=jpoint(nbr)
          ip1=min(i1+1,jmax1)
          im1=max(i1-1,1)
          jp1=min(j1+1,kmax1)
          jm1=max(j1-1,1)

          x_1(nbr,id)  = x1(i1,j1)
          y_1(nbr,id)  = y1(i1,j1)
          rhoi         = 1.0/q1(i1,j1,1)
          Q_1(1,nbr,id) = q1(i1,j1,1)*q1(i1,j1,5)
          Q_1(2,nbr,id) = q1(i1,j1,2)*rhoi
          Q_1(3,nbr,id) = q1(i1,j1,3)*rhoi
          Q_1(4,nbr,id) = gm1*( q1(i1,j1,4)-0.5*(q1(i1,j1,2)**2+q1(i1,j1,3)**2)
     <                                         /q1(i1,j1,1))*q1(i1,j1,5)

c         Gradients

          rhoi          = 1.0/q1(ip1,j1,1)
          qip1 		= q1(ip1,j1,1)*q1(ip1,j1,5)
          qip2          = q1(ip1,j1,2)*rhoi
          qip3		= q1(ip1,j1,3)*rhoi
          qip4          = gm1*( q1(ip1,j1,4)-0.5*(q1(ip1,j1,2)**2+q1(ip1,j1,3)**2)
     <                                         /q1(ip1,j1,1))*q1(ip1,j1,5)
          rhoi          = 1.0/q1(im1,j1,1)
          qim1 		= q1(im1,j1,1)*q1(im1,j1,5)
          qim2          = q1(im1,j1,2)*rhoi
          qim3		= q1(im1,j1,3)*rhoi
          qim4          = gm1*( q1(im1,j1,4)-0.5*(q1(im1,j1,2)**2+q1(im1,j1,3)**2)
     <                                         /q1(im1,j1,1))*q1(im1,j1,5)

          rhoi          = 1.0/q1(i1,jp1,1)
          qjp1 		= q1(i1,jp1,1)*q1(i1,jp1,5)
          qjp2          = q1(i1,jp1,2)*rhoi
          qjp3		= q1(i1,jp1,3)*rhoi
          qjp4          = gm1*( q1(i1,jp1,4)-0.5*(q1(i1,jp1,2)**2+q1(i1,jp1,3)**2)
     <                                         /q1(i1,jp1,1))*q1(i1,jp1,5)
          rhoi          = 1.0/q1(i1,jm1,1)
          qjm1 		= q1(i1,jm1,1)*q1(i1,jm1,5)
          qjm2          = q1(i1,jm1,2)*rhoi
          qjm3		= q1(i1,jm1,3)*rhoi
          qjm4          = gm1*( q1(i1,jm1,4)-0.5*(q1(i1,jm1,2)**2+q1(i1,jm1,3)**2)
     <                                         /q1(i1,jm1,1))*q1(i1,jm1,5)

          Qx_1(1,nbr,id)=0.5*((qip1-qim1)*xx1(i1,j1)+(qjp1-qjm1)*yx1(i1,j1))
          Qx_1(2,nbr,id)=0.5*((qip2-qim2)*xx1(i1,j1)+(qjp2-qjm2)*yx1(i1,j1))
          Qx_1(3,nbr,id)=0.5*((qip3-qim3)*xx1(i1,j1)+(qjp3-qjm3)*yx1(i1,j1))
          Qx_1(4,nbr,id)=0.5*((qip4-qim4)*xx1(i1,j1)+(qjp4-qjm4)*yx1(i1,j1))

          Qy_1(1,nbr,id)=0.5*((qip1-qim1)*xy1(i1,j1)+(qjp1-qjm1)*yy1(i1,j1))
          Qy_1(2,nbr,id)=0.5*((qip2-qim2)*xy1(i1,j1)+(qjp2-qjm2)*yy1(i1,j1))
          Qy_1(3,nbr,id)=0.5*((qip3-qim3)*xy1(i1,j1)+(qjp3-qjm3)*yy1(i1,j1))
          Qy_1(4,nbr,id)=0.5*((qip4-qim4)*xy1(i1,j1)+(qjp4-qjm4)*yy1(i1,j1))

          enddo

      enddo

c        Send stuff to Praveendora's Magic Box

         call meshless(dt,ndonor,x_0,y_0,Q_0,x_1,y_1,Q_1,Qx_1,Qy_1)

c        Update actually happens here. Q_0 is just dt*Residue

      do id=1,ndonor
         
         i=imesh(id,1)
         j=imesh(id,2)

         q(i,j,1)=q(i,j,1)+Q_0(1,id)/q(i,j,5)
         q(i,j,2)=q(i,j,2)+Q_0(2,id)/q(i,j,5)
         q(i,j,3)=q(i,j,3)+Q_0(3,id)/q(i,j,5)
         q(i,j,4)=q(i,j,4)+Q_0(4,id)/q(i,j,5)

      enddo

      return
      end





