#######################################################################
#set size square
set auto
set nokey
set contour
set param
set nosurface
set view 0,360
#set cntrparam levels 25
set cntrparam levels incremental 0, 0.05, 1
set term table
set out 'cnt.dat'
set xran[0:1]
set yran[0:1]
set size ratio -1
splot 'result.dat' u 1:2:3 w l
#set term x11
set term gif
set out 'cont0.gif'
#set noxtics
#set noytics
p 'cnt.dat' w l
########################################################################
