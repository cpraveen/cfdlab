# Linear advection equation in 1d

Solves linear advection eqn using upwind finite volume scheme, weno5 reconstruction and periodic bc.

We use PETSc vectors to manage parallelization but there is no matrix that is being solved.

To run:
```
make
mpirun -np 4 ./convect
```
Then plot solution in gnuplot
```
gnuplot> p "sol0.dat" w l,"sol.dat" w p
```
or run the python script plt.py

By default, the code uses 200 grid points. To specify number of grid points,
```
mpirun -np 4 ./convect -da_grid_x 1000
```
