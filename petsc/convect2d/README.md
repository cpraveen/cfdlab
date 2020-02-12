# 2d linear advection equation

This code solves
```
u_t + u_x + u_y = 0
```
with periodic boundary conditions. To run this code, do
```
make
rm -f sol*.plt
mpirun -np 4 ./convect
sh ./merge.sh
```
Open the `sol*.plt` files in VisIt and animate the solution. Note that it is important to delete the old solution plt files.

You can specify the grid size as a command line argument, below we specify a 100x100 mesh
```
rm -f sol*.plt
mpirun -np 4 ./convect -da_grid_x 100 -da_grid_y 100
sh ./merge.sh
```
