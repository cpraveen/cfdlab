# 2d Euler solution using central finite volume

This code solves 2d Euler equations on Cartesian mesh using central finite volume method with periodic boundary conditions.

This makes use of time stepping schemes in Petsc. To solve du/dt = R(t,u) you must implement R inside the function RHSFunction. Specify either dt or cfl. If both are given, then cfl will be used to compute time step. 

Some options
```
-problem  vortex | density
-flux     central | kepec | kep | mkep | kg | ducros
-order    2 | 4
```

The code requires PETSc to compile. Set the variable
```
export PETSC_DIR=/path/to/your/petsc/installation
```
in your shell.

First compile the code
```
make
```
If you dont specify any scheme
```
rm -f sol*.plt
mpirun -np 4 ./ts -problem vortex -flux kepec -order 2 \
                  -da_grid_x 100 -da_grid_y 100 \
                  -Tf 20.0 -cfl 1.8 -si 100 -ts_monitor 
sh ./merge.sh
```
it will use 2-stage, 2-nd order SSPRK scheme.

You can open the plt files using Tecplot or VisIt.

The following will use 4-stage, 3-order SSPRK scheme.
```
rm -f sol*.plt
mpirun -np 4 ./ts -problem vortex -flux kepec -order 2 \
                  -da_grid_x 100 -da_grid_y 100 \
                  -Tf 20.0 -cfl 1.8 -si 100 \
                  -ts_type ssp -ts_ssp_type rks3 -ts_ssp_nstages 4 -ts_monitor 
sh ./merge.sh
```
To use the classical RK4 scheme
```
rm -f sol*.plt
mpirun -np 4 ./ts -problem vortex -flux kepec -order 2 \
                  -da_grid_x 100 -da_grid_y 100 \
                  -Tf 20.0 -cfl 0.8 -si 100 \
                  -ts_type rk -ts_rk_type 4 -ts_adapt_type none -ts_monitor 
sh ./merge.sh
```