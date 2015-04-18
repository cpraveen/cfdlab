set xlabel 'x'
set ylabel 'f_1, f_2'
p 'sol.dat' u 1:2 t 'f1' w l, \
  'sol.dat' u 1:3 t 'f2' w l
