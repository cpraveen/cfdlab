set xlabel 'x'
set ylabel 'f_1, f_2'
p 'sol.dat' u 1:2 t 'f1' w l lt 2, \
  'sol.dat' u 1:3 t 'f2' w l lt 3
set term postscript enhanced color
set out 'sol.eps'
replot
set term x11
set out
