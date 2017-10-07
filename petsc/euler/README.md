# 2d Euler solution using WENO5 finite volume

This code solves 2d Euler equations on Cartesian mesh using WENO5 finite volume method with periodic boundary conditions.
```
make ssprk
rm -f sol*.plt
mpirun -np 4 ./ssprk -da_grid_x 100 -da_grid_y 100 -Tf 10.0 -cfl 0.8 -si 100
sh ./merge.sh
```
You can open the plt files using Tecplot of VisIt.

## TS version

This makes use of time stepping schemes in Petsc. To solve du/dt = R(t,u) you must implement R inside the function RHSFunction. Specify either dt or cfl. If both are given, then cfl will be used to compute time step. The following will use 4-stage, 3-order SSPRK scheme.

```
make ts
rm -f sol*.plt
mpirun -np 4 ./ts -da_grid_x 100 -da_grid_y 100 -Tf 10.0 -cfl 1.8 -si 100 \
                  -ts_type ssp -ts_ssp_type rks3 -ts_ssp_nstages 4 -ts_monitor 
sh ./merge.sh
```
