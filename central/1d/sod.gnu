set term postscript enhanced
set out 'sod.eps'

set xlabel 'x'
set ylabel 'Density'
p 'sol.txt' u 1:2 t 'KT' w lp,'sod.txt' u 1:2 t 'Exact' w l

set xlabel 'x'
set ylabel 'Velocity'
p 'sol.txt' u 1:3 t 'KT' w lp,'sod.txt' u 1:3 t 'Exact' w l

set xlabel 'x'
set ylabel 'Pressure'
p 'sol.txt' u 1:4 t 'KT' w lp,'sod.txt' u 1:4 t 'Exact' w l
