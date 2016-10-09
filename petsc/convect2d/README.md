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
Open the `sol*.plt` files in VisIt and animate the solution.
