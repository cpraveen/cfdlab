# 2d Euler solution using WENO5 finite volume

## SSPRK version (ssprk.c, finite volume WENO)
This code solves 2d Euler equations on Cartesian mesh for the isentropic vortex problem using WENO5 finite volume method with periodic boundary conditions.

```shell
make ssprk
rm -f sol*.plt
mpirun -np 4 ./ssprk -da_grid_x 100 -da_grid_y 100 -Tf 20.0 -cfl 0.8 -si 100
sh ./merge.sh
visit -o sol*.plt
```

You can open the plt files using Tecplot of VisIt.

## TS version (ts.c, finite volume WENO)

This is similar to `ssprk.c` but makes use of time stepping schemes in Petsc. To solve du/dt = R(t,u) you must implement R inside the function RHSFunction. Specify either dt or cfl. If both are given, then cfl will be used to compute time step. First compile the code

```shell
make ts
```

If you dont specify any scheme

```shell
rm -f sol*.plt
mpirun -np 4 ./ts -da_grid_x 100 -da_grid_y 100 -Tf 20.0 -cfl 1.8 -si 100 \
                  -ts_monitor
sh ./merge.sh
visit -o sol*.plt
```

it will use 2-stage, 2-nd order SSPRK scheme.

The following will use 4-stage, 3-order SSPRK scheme.

```shell
rm -f sol*.plt
mpirun -np 4 ./ts -da_grid_x 100 -da_grid_y 100 -Tf 20.0 -cfl 1.8 -si 100 \
                  -ts_type ssp -ts_ssp_type rks3 -ts_ssp_nstages 4 -ts_monitor
sh ./merge.sh
visit -o sol*.plt
```

To use the classical RK4 scheme

```shell
rm -f sol*.plt
mpirun -np 4 ./ts -da_grid_x 100 -da_grid_y 100 -Tf 20.0 -cfl 0.8 -si 100 \
                  -ts_type rk -ts_rk_type 4 -ts_adapt_type none -ts_monitor
sh ./merge.sh
visit -o sol*.plt
```

## TS version (fdweno.c, finite difference WENO)

This is the most sophisticated of the three codes.  Compile the code

```shell
make fdweno PROBLEM=ISENTROPIC WENO=z
```

Run this similar to ts.c code. Other options for PROBLEM are

* ISENTROPIC: Isentropic vortex with periodic bc
* SHOCKREF  : Shock reflection (has steady solution for large time)
* SHOCKVORTEX: Shock vortex interaction
* RIEMANN2D: 2-D Riemann problem
* KH: Kelvin-Helmholtz instability

Options for WENO are js and z. You can see some make options

```shell
make help
```

### Isentropic vortex

```shell
rm -f fdweno && make fdweno PROBLEM=ISENTROPIC WENO=z
rm -rf sol*.plt
mpirun -np 4 ./fdweno -da_grid_x 100 -da_grid_y 100 -Tf 5.0 -cfl 0.8 -si 100 \
       -ts_type ssp -ts_ssp_type rks3 -ts_ssp_nstages 4 -ts_monitor
sh ./merge.sh
visit -o sol*.plt
```

### Shock reflection

```shell
rm -f fdweno && make fdweno PROBLEM=SHOCKREF WENO=z
rm -rf sol*.plt
mpirun -np 4 ./fdweno -da_grid_x 200 -da_grid_y 50 -Tf 5.0 -cfl 0.8 -si 100 \
       -ts_type ssp -ts_ssp_type rks3 -ts_ssp_nstages 4 -ts_monitor
sh ./merge.sh
visit -o sol*.plt
```

### Shock vortex interaction

```shell
rm -f fdweno && make fdweno PROBLEM=SHOCKVORTEX WENO=z
rm -rf sol*.plt
mpirun -np 4 ./fdweno -da_grid_x 100 -da_grid_y 100 -Tf 0.5 -cfl 0.8 -si 10 \
       -ts_type ssp -ts_ssp_type rks3 -ts_ssp_nstages 4 -ts_monitor
sh ./merge.sh
visit -o sol*.plt
```

### Double mach reflection

```shell
rm -f fdweno && make fdweno PROBLEM=DMR WENO=z
rm -rf sol*.plt
mpirun -np 4 ./fdweno -da_grid_x 400 -da_grid_y 100 -Tf 0.2 -cfl 0.8 -si 100 \
       -ts_type ssp -ts_ssp_type rks3 -ts_ssp_nstages 4 -ts_monitor
sh ./merge.sh
visit -o sol*.plt
```

### 2-D Riemann problem

```shell
rm -f fdweno && make fdweno PROBLEM=2DRIEMANN WENO=z
rm -rf sol*.plt
mpirun -np 4 ./fdweno -da_grid_x 100 -da_grid_y 100 -Tf 0.8 -cfl 0.8 -si 100 \
       -ts_type ssp -ts_ssp_type rks3 -ts_ssp_nstages 4 -ts_monitor
sh ./merge.sh
visit -o sol*.plt
```
