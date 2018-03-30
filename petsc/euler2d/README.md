# 2d Euler solution using WENO5 finite volume

This code solves 2d Euler equations on Cartesian mesh using WENO5 finite volume method with periodic boundary conditions.
```
make ssprk
rm -f sol*.plt
mpirun -np 4 ./ssprk -da_grid_x 100 -da_grid_y 100 -Tf 20.0 -cfl 0.8 -si 100
sh ./merge.sh
```
You can open the plt files using Tecplot of VisIt.

## TS version (ts.c, finite volume WENO)

This makes use of time stepping schemes in Petsc. To solve du/dt = R(t,u) you must implement R inside the function RHSFunction. Specify either dt or cfl. If both are given, then cfl will be used to compute time step. First compile the code
```
make ts
```
If you dont specify any scheme
```
rm -f sol*.plt
mpirun -np 4 ./ts -da_grid_x 100 -da_grid_y 100 -Tf 20.0 -cfl 1.8 -si 100 \
                  -ts_monitor 
sh ./merge.sh
```
it will use 2-stage, 2-nd order SSPRK scheme.

The following will use 4-stage, 3-order SSPRK scheme.
```
rm -f sol*.plt
mpirun -np 4 ./ts -da_grid_x 100 -da_grid_y 100 -Tf 20.0 -cfl 1.8 -si 100 \
                  -ts_type ssp -ts_ssp_type rks3 -ts_ssp_nstages 4 -ts_monitor 
sh ./merge.sh
```
To use the classical RK4 scheme
```
rm -f sol*.plt
mpirun -np 4 ./ts -da_grid_x 100 -da_grid_y 100 -Tf 20.0 -cfl 0.8 -si 100 \
                  -ts_type rk -ts_rk_type 4 -ts_adapt_type none -ts_monitor 
sh ./merge.sh
```

## TS version (fdweno.c, finite difference WENO)

This is the most sophisticated of the three codes.  Compile the code
```
make fdweno PROBLEM=ISENTROPIC
```
Run this similar to ts.c code. Other options for PROBLEM are

 * ISENTROPIC: Isentropic vortex with periodic bc
 * SHOCKREF  : Shock reflection
 * SHOCKVORTEX: Shock vortex interaction
