= Linear advection equation in 2d using dmplex

Make grid
```
gmsh -2 -format msh2 tri.geo
```

Compile code
```
make
```

Run it
```
./convect -ts_monitor -ts_view
```
