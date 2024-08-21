reset

set term qt 1
set title "Density"
p "sol_con.txt" u 1:2 w l,\
  "sol_noncon1.txt" u 1:2 w l, \
  "sol_noncon2.txt" u 1:2 w l
set xlabel "x"
set ylabel "Density"

set term qt 2
set title "Velocity"
p "sol_con.txt" u 1:3 w l,\
  "sol_noncon1.txt" u 1:3 w l, \
  "sol_noncon2.txt" u 1:3 w l
set xlabel "x"
set ylabel "Velocity"

set term qt 3
set title "Pressure"
p "sol_con.txt" u 1:4 w l,\
  "sol_noncon1.txt" u 1:4 w l, \
  "sol_noncon2.txt" u 1:4 w l
set xlabel "x"
set ylabel "Pressure"
