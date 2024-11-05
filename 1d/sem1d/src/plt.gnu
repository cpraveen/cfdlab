set term pdf
set out 'sol.pdf'
plot 'sol0.dat' t 'Initial condition' w lp lw 2 pt 2,\
     'sol.dat'  t 'SEM' w lp lw 2 pt 3,\
     'ex.dat'   t 'Exact' w lp lw 2 pt 4
