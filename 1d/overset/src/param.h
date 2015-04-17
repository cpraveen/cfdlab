      integer ngridmax, npmax, nolapmax
      parameter(ngridmax=10, npmax=1000, nolapmax=1000)
      integer FIXED, UPDATE, OVRLAP
      parameter(FIXED=0, UPDATE=1, OVRLAP=2)

      integer ngrid, nopts
      integer np(ngridmax), utype(ngridmax,npmax)
      integer olap(7,nolapmax)
      real    x(ngridmax,npmax), xmin(ngridmax), xmax(ngridmax),
     &        dx(ngridmax)

      common/ivars/ngrid,nopts,np,utype,olap
      common/rvars/x,xmin,xmax,dx
