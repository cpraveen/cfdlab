reset 
set term postscript color
set out 'result.eps'
set grid

p 'shape.dat' u 1:2 t 'Nozzle shape' w l lw 2

p 'flow.dat' u 1:2 t 'Density' w lp pt 6

p  'flow.dat' u 1:3 t 'Velocity' w lp pt 6

p  'flow.dat' u 1:4 t 'Pressure' w lp pt 6

p  'flow.dat' u 1:5 t 'Mach' w lp pt 6

p 'flowb.dat' u 1:2 t 'Adjoint density' w lp pt 6

p  'flowb.dat' u 1:3 t 'Adjoint momentum' w lp pt 6

p  'flowb.dat' u 1:4 t 'Adjoint energy' w lp pt 6
