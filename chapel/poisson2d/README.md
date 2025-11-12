# Solve 2d Poisson equation

Uses V-cycle multigrid. It requires the VTK writer, see the makefile to set the path.

Compile the code

```shell
make
```

Run the code

```shell
./main --nx 128 --ny 128 --levels 7
python ./plot.py
```

Increase grid size

```shell
./main --nx 256 --ny 256 --levels 8
python ./plot.py
```

For best convergence rate, choose `nx = ny = 2^levels`.
