      subroutine ReadParams
      implicit none
      include 'vars.h'

      integer fid

      fid = 50

      open(fid, file='param.in')
      read(fid,*)xmin, xmax, ymin, ymax
      read(fid,*)nx, ny
      read(fid,*)perturb, pertx, perty
      read(fid,*)order
      read(fid,*)kkk
      read(fid,*)ictype
      read(fid,*)veltype
      read(fid,*)dt
      read(fid,*)Ttime
      read(fid,*)makeanim, animinterval

      return
      end
