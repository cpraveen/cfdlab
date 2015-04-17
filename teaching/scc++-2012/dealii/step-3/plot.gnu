set hidden3d
splot "solution.gpl" w l

pause 2

set hidden3d
set pm3d at b
replot

pause 2

set pm3d map
set size ratio -1
replot
