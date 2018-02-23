out="error.txt"

# Grid sizes
N="50 100 200 400"

rm -f $out
touch $out

for n in $N
do
   echo "n =" $n
   mpirun -np 4 ./fdweno -da_grid_x $n -da_grid_y $n -Tf 20.0 -cfl 0.8 \
          -si 100000 -ts_type rk -ts_rk_type 4 -ts_adapt_type none \
          > log.txt
   tail -n1 log.txt
   tail -n1 log.txt >> $out
done
