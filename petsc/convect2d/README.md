# 2d linear advection equation

This code solves

```shell
u_t + u_x + u_y = 0
```

with periodic boundary conditions. To run this code, do

```shell
make
rm -f sol*.plt
mpirun -np 4 ./convect
sh ./merge.sh
```

Open the `sol*.plt` files in VisIt and animate the solution.

```shell
visit -o sol*.plt
```

Note that it is important to delete the old solution plt files before running the code.

You can specify the grid size as a command line argument, below we specify a 100x100 mesh

```shell
rm -f sol*.plt
mpirun -np 4 ./convect -da_grid_x 100 -da_grid_y 100
sh ./merge.sh
visit -o sol*.plt
```
