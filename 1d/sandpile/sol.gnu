set term pdf
set out 'sol.pdf'

set grid
set ylabel 'u'
set xlabel 'x'
p 'sol.dat' u 1:2 t 'DG' w l lw 2,\
  'exact.dat' u 1:2 t 'Exact' w l lw 2 lc 3

set grid
set ylabel 'v'
set xlabel 'x'
p 'sol.dat' u 1:3 t 'DG' w l lw 2,\
  'exact.dat' u 1:3 t 'Exact' w l lw 2 lc 3
