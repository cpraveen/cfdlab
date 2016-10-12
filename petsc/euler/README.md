# 2d Euler solution using WENO5 finite volume

This code solves 2d Euler equations on Cartesian mesh using WENO5 finite volume method with periodic boundary conditions.
```
make
rm -f sol*.plt
mpirun -np 4 ./euler -da_grid_x 100 -da_grid_y 100 -Tf 10.0 -cfl 0.4
sh ./merge.sh
```
You can open the plt files using Tecplot of VisIt.
